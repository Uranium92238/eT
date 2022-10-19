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
module cnto_tool_class
!
!!
!!    Correlated natural transition orbital tool class
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Handles construction and diagonalization of M and N
!!
!!    M_ij += sum_a X_ai X_aj
!!    M_ij += 1/2 sum_abl(1 + delta_ai,bl delta_i,j) X_aibl X_ajbl)
!!
!!    N_ab += sum_i X_ai X_bi
!!    N_ab += 1/2 sum_cij(1 + delta_ai,cj delta_a,b) X_aicj X_bicj)
!!
!
   use parameters
   use memory_manager_class, only: mem
   use nto_tool_class,       only: nto_tool
!
   implicit none
!
   type, extends(nto_tool) :: cnto_tool
!
   contains
!
      procedure, public :: add_contributions_to_M_and_N &
                        => add_contributions_to_M_and_N_cnto_tool
!
      procedure, public :: add_doubles_to_M_and_N &
                        => add_doubles_to_M_and_N_cnto_tool
!
   end type cnto_tool
!
   interface cnto_tool
!
      procedure :: new_cnto_tool
!
   end interface cnto_tool
!
contains
!
!
   function new_cnto_tool(n_o, n_v, n_ao, X_length) result(this)
!!
!!    New cnto tool
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      type(cnto_tool) :: this
!
      integer, intent(in) :: n_o, n_v, n_ao, X_length
!
      this%n_ao = n_ao
      this%n_o = n_o
      this%n_v = n_v
      this%X_length = X_length
!
      call this%prepare('cnto')
!
   end function new_cnto_tool
!
!
   subroutine add_doubles_to_M_and_N_cnto_tool(this, X2)
!!
!!    Add doubles to M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Add the doubles contribution of an excited state X to M and N
!!
      implicit none
!
      class(cnto_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_v, this%n_o, this%n_v, this%n_o), intent(in) :: X2
!
      integer :: a, i
!
      call dgemm('T', 'N',             &
                 this%n_o,             &
                 this%n_o,             &
                 this%n_v**2*this%n_o, &
                 half,                 &
                 X2,                   & ! R_bla_i
                 this%n_v**2*this%n_o, &
                 X2,                   & ! R_bla_j
                 this%n_v**2*this%n_o, &
                 one,                  &
                 this%M,               &
                 this%n_o)
!
!$omp parallel do private(a, i)
      do i = 1, this%n_o
         do a = 1, this%n_v
!
            this%M(i, i) = this%M(i, i) + half*X2(a,i,a,i)**2
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'T',             &
                 this%n_v,             &
                 this%n_v,             &
                 this%n_v*this%n_o**2, &
                 half,                 &
                 X2,                   &  ! R_a_icj
                 this%n_v,             &
                 X2,                   &  ! R_b_icj
                 this%n_v,             &
                 one,                  &
                 this%N,               &
                 this%n_v)
!
!$omp parallel do private(a, i)
      do a = 1, this%n_v
         do i = 1, this%n_o
!
            this%N(a, a) = this%N(a, a) + half*X2(a,i,a,i)**2
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_doubles_to_M_and_N_cnto_tool
!
!
   subroutine add_contributions_to_M_and_N_cnto_tool(this, X)
!!
!!    Add contributions to M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!
      use reordering, only: squareup
!
      implicit none
!
      class(cnto_tool), intent(inout) :: this
!
      real(dp), dimension(this%X_length), intent(in) :: X
      real(dp), dimension(:,:,:,:), allocatable :: X2
!
      call this%add_singles_to_M_and_N(X(1 : this%n_o * this%n_v))
!
      call mem%alloc(X2, this%n_v, this%n_o, this%n_v, this%n_o)
!
      call squareup(X(this%n_o * this%n_v + 1 :), X2, this%n_o * this%n_v)
!
      call this%add_doubles_to_M_and_N(X2)
!
      call mem%dealloc(X2, this%n_v, this%n_o, this%n_v, this%n_o)
!
   end subroutine add_contributions_to_M_and_N_cnto_tool
!
!
end module cnto_tool_class
