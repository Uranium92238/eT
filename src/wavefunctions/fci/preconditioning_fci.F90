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
submodule (fci_class) preconditioning_fci
!
!!
!! FCI preconditioning submodule
!!
!! Gathers routines to precondition a FCI calculation
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_h_diagonal(wf, h_diagonal)
!!
!!    Construct h diagonal
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci), intent(inout) :: wf
      real(dp), dimension(wf%n_alpha_strings, wf%n_beta_strings), intent(out) :: h_diagonal
!
      integer :: p, q, i, j, alpha_string, beta_string
!
      real(dp) :: h_one_electron, h_two_electron, h_nuclear_repulsion
!
      integer, dimension(:,:), allocatable :: occupied_alpha, occupied_beta
!
      type(timings), allocatable :: timer
!
      timer = timings('Construct Hamiltonian diagonal', pl='n')
      call timer%turn_on
!
      h_nuclear_repulsion = wf%get_nuclear_repulsion()
!
      call mem%alloc(occupied_alpha, wf%n_alpha_strings, wf%n_alpha)
      call mem%alloc(occupied_beta,  wf%n_beta_strings,  wf%n_beta)
!
      call wf%calculate_occupied_indices_in_strings(wf%n_alpha_strings, &
                                                    wf%n_mo,            &
                                                    wf%n_alpha,         &
                                                    wf%alpha_strings,   &
                                                    occupied_alpha)
!
      call wf%calculate_occupied_indices_in_strings(wf%n_beta_strings,  &
                                                    wf%n_mo,            &
                                                    wf%n_beta,          &
                                                    wf%beta_strings,    &
                                                    occupied_beta)
!
!$omp parallel do private(alpha_string, beta_string, i, j, h_one_electron, h_two_electron)
      do alpha_string = 1, wf%n_alpha_strings
         do beta_string = 1, wf%n_beta_strings
!
            h_one_electron = zero
            h_two_electron = zero
!
            do p = 1, wf%n_alpha
!
               i = occupied_alpha(alpha_string, p)
               h_one_electron = h_one_electron + wf%h_pq(i,i)
!
               do q = 1, wf%n_alpha
!
                  j = occupied_alpha(alpha_string, q)
                  h_two_electron = h_two_electron &
                                 + (wf%g_pqrs(i,i,j,j) &
                                 -  wf%g_pqrs(i,j,j,i))
!
               end do
!
               do q = 1, wf%n_beta
!
                  j = occupied_beta(beta_string, q)
                  h_two_electron = h_two_electron &
                                 + two*wf%g_pqrs(i,i,j,j)
!
               end do
!
            enddo
!
            do p = 1, wf%n_beta
!
               i = occupied_beta(beta_string, p)
               h_one_electron = h_one_electron + wf%h_pq(i,i)
!
               do q = 1, wf%n_beta
!
                  j = occupied_beta(beta_string, q)
                  h_two_electron = h_two_electron &
                                 + (wf%g_pqrs(i,i,j,j) &
                                 -  wf%g_pqrs(i,j,j,i))
!
               end do
!
            enddo
!
            h_diagonal(alpha_string, beta_string) = h_one_electron &
                                                  + half*h_two_electron &
                                                  + h_nuclear_repulsion
!
         end do
      end do
!$omp end parallel do
!
      call mem%dealloc(occupied_alpha, wf%n_alpha_strings, wf%n_alpha)
      call mem%dealloc(occupied_beta,  wf%n_beta_strings,  wf%n_beta)
!
      call timer%turn_off
!
   end subroutine construct_h_diagonal
!
!
end submodule preconditioning_fci
