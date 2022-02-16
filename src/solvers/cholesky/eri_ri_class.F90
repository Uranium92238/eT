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
module eri_ri_class
!
!!
!!    Electron repulsion integral (ERI) resolution-of-identity (RI) class
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!!    Handles the generation of RI vectors (in the MO basis)
!!
!!       L_pq^J = sum_{K} [Q^-1]_JK (K|pq)
!!
!!    where J, K are elements of the auxiliary basis and
!!    Q is the Cholesky factor of the matrix S_JK = (J|K).
!!
!
   use parameters
!
   use ao_tool_class, only : ao_tool
   use ri_basis_class, only : ri_basis
   use memory_manager_class, only: mem
   use global_out, only: output
!
   implicit none
!
   type :: eri_ri
!
      real(dp), dimension(:,:), allocatable, private :: Q_inverse, Q, S
!
      type(ri_basis),  private :: aux
!
   contains
!
      procedure, public :: initialize &
                        => initialize_eri_ri
!
      procedure, public :: run &
                        => run_eri_ri
!
      procedure, public :: construct_cholesky_mo_vectors &
                        => construct_cholesky_mo_vectors_eri_ri

!
      procedure, public :: get_n_J &
                        => get_n_J_eri_ri

      procedure, public :: cleanup &
                        => cleanup_eri_ri
!
      procedure, private :: invert_Q
      procedure, private :: construct_S
      procedure, private :: decompose_S
      procedure, private :: print_summary
!
   end type eri_ri
!
!
   interface eri_ri
!
      procedure :: new_eri_ri
!
   end interface eri_ri
!
contains
!
   function new_eri_ri(basis_set, eri_precision) result(this)
!!
!!    New eri ri
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!
      use iso_c_binding, only: c_char, c_null_char, c_double, c_int
!
      implicit none
!
      character(len=200), intent(in) :: basis_set
      real(dp), intent(in) :: eri_precision
!
      type(eri_ri) :: this
!
      this%aux = ri_basis(basis_set, eri_precision)
!
   end function new_eri_ri
!
!
   subroutine initialize_eri_ri(this)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      class(eri_ri), intent(inout) :: this
!
      call this%aux%initialize_shells()
!
      call mem%alloc(this%S, this%aux%n_ao, this%aux%n_ao)
      call mem%alloc(this%Q, this%aux%n_ao, this%aux%n_ao)
      call mem%alloc(this%Q_inverse, this%aux%n_ao, this%aux%n_ao)
!
   end subroutine initialize_eri_ri
!
!
   subroutine run_eri_ri(this, ao)
!!
!!    Run
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      class(eri_ri), intent(inout) :: this
      type(ao_tool), intent(in) :: ao
!
      call this%construct_S(ao)
      call this%decompose_S()
      call this%invert_Q()
!
      call this%print_summary()
!
   end subroutine run_eri_ri
!
!
   subroutine construct_cholesky_mo_vectors_eri_ri(this, ao, n_ao, n_mo, &
                                                   orbital_coefficients, cd_tool)
!!
!!    Construct Cholesky MO vectors
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!!    Constructs vectors
!!
!!       L_pq^J = sum_K [Q^-1]_JK (K|pq)
!!              = sum_K sum_wx [Q^-1]_JK (K|wx) C_wp C_xq
!!
!!    these are not Cholesky vectors since the auxilliary basis is
!!    the RI basis.
!!
      use array_utilities, only: zero_array
      use reordering, only: sort_123_to_132
      use batching_index_class, only: batching_index
      use abstract_eri_cholesky_class, only: abstract_eri_cholesky
!
      implicit none
!
      class(eri_ri) :: this
!
      type(ao_tool) :: ao
!
      integer, intent(in) :: n_mo, n_ao
!
      real(dp), dimension(n_ao, n_mo), intent(in) :: orbital_coefficients
!
      class(abstract_eri_cholesky), intent(inout) :: cd_tool
!
      real(dp), dimension(:,:,:), allocatable :: X_Kwq, X_Kqw, X_Kqp, L_Jqp, g_Kwx, X_Kqp_red, L_Jpq
      real(dp), dimension(ao%max_sh_size**2*this%aux%max_dim), target :: g_KCD
      real(dp), dimension(:,:,:), pointer :: g_KCD_p
!
      integer :: J, C, D, w, x, z, i
      integer :: req0, req1, current_q_batch

      type(batching_index) :: batch_q
!
      req0 = this%aux%max_dim*n_ao**2
      req1 = 2*n_mo*this%aux%n_ao + max(n_ao*this%aux%max_dim*2, n_mo*this%aux%max_dim*2)
!
      batch_q = batching_index(n_mo)
!
      call mem%batch_setup(batch_q, req0, req1, tag='RI Cholesky vectors')
!
      do current_q_batch = 1, batch_q%num_batches
!
         call batch_q%determine_limits(current_q_batch)
!
         call mem%alloc(X_Kqp, this%aux%n_ao, batch_q%length, n_mo)
!
         do J = 1, this%aux%n_shells
!
            call mem%alloc(g_Kwx, this%aux%get_shell_length(J), n_ao, n_ao)
!
!$omp parallel do private(C, D, w, x, z, g_KCD_p, g_KCD)
            do C = 1, ao%n_sh
               do D = 1, ao%n_sh
!
               call ao%get_eri_3c(g_KCD, this%aux%get_shell(J), C, D)
!
               g_KCD_p(1 : this%aux%get_shell_length(J), &
                     1 : ao%shells(C)%length, &
                     1 : ao%shells(D)%length) &
                     => g_KCD(1 : &
                           this%aux%get_shell_length(J)*ao%shells(C)%length*ao%shells(D)%length)
!
               do x = ao%shells(D)%first, ao%shells(D)%get_last()
                  do w = ao%shells(C)%first, ao%shells(C)%get_last()
                     do z = 1, this%aux%get_shell_length(J)
!

                        g_Kwx(z, w, x) = g_KCD_p(z, &
                                             w - ao%shells(C)%first + 1, &
                                             x - ao%shells(D)%first + 1)
!
                     enddo
                  enddo
               enddo
            enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(X_Kwq, this%aux%get_shell_length(J), n_ao, batch_q%length)
!
            call dgemm('N', 'N',                                 &
                        n_ao*this%aux%get_shell_length(J),       &
                        batch_q%length,                          &
                        n_ao,                                    &
                        one,                                     &
                        g_Kwx,                                   &
                        n_ao*this%aux%get_shell_length(J),       &
                        orbital_coefficients(1, batch_q%first),  &
                        n_ao,                                    &
                        zero,                                    &
                        X_Kwq,                                   &
                        n_ao*this%aux%get_shell_length(J))
!
            call mem%dealloc(g_Kwx, this%aux%get_shell_length(J), n_ao, n_ao)
!
            call mem%alloc(X_Kqw, this%aux%get_shell_length(J), batch_q%length, n_ao)
            call sort_123_to_132(X_Kwq, X_Kqw, this%aux%get_shell_length(J), n_ao, batch_q%length)
            call mem%dealloc(X_Kwq, this%aux%get_shell_length(J), n_ao, batch_q%length)
!
            call mem%alloc(X_Kqp_red, this%aux%get_shell_length(J), batch_q%length, n_mo)
!
            call dgemm('N', 'N',                                       &
                        batch_q%length*this%aux%get_shell_length(J),   &
                        n_mo,                                          &
                        n_ao,                                          &
                        one,                                           &
                        X_Kqw,                                         &
                        batch_q%length*this%aux%get_shell_length(J),   &
                        orbital_coefficients,                          &
                        n_ao,                                          &
                        zero,                                          &
                        X_Kqp_red,                                     &
                        batch_q%length*this%aux%get_shell_length(J))
!
            call mem%dealloc(X_Kqw, this%aux%get_shell_length(J), batch_q%length, n_ao)
!
            do i = this%aux%get_shell_first(J), this%aux%get_shell_last(J)
!
               X_Kqp(i, :, :) = X_Kqp_red(i - this%aux%get_shell_first(J) + 1, :, :)
!
            enddo
!
            call mem%dealloc(X_Kqp_red, this%aux%get_shell_length(J), batch_q%length, n_mo)
!
         enddo
!
         call mem%alloc(L_Jqp, this%aux%n_ao, batch_q%length, n_mo)
!
         call dgemm('N', 'N',             &
                     this%aux%n_ao,       &
                     n_mo*batch_q%length, &
                     this%aux%n_ao,       &
                     one,                 &
                     this%Q_inverse,      &
                     this%aux%n_ao,       &
                     X_Kqp,               &
                     this%aux%n_ao,       &
                     zero,                &
                     L_Jqp,               &
                     this%aux%n_ao)
!
         call mem%dealloc(X_Kqp, this%aux%n_ao, batch_q%length, n_mo)
         call mem%alloc(L_Jpq, this%aux%n_ao, n_mo, batch_q%length)
         call sort_123_to_132(L_Jqp, L_Jpq, this%aux%n_ao, batch_q%length, n_mo)
         call mem%dealloc(L_Jqp, this%aux%n_ao, batch_q%length, n_mo)
!
         call cd_tool%set(L_Jpq, 1, n_mo, batch_q%first, batch_q%get_last())
!
         call mem%dealloc(L_Jpq, this%aux%n_ao, n_mo, batch_q%length)
!
      enddo
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_mo_vectors_eri_ri
!
!
   subroutine cleanup_eri_ri(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      class(eri_ri) :: this
!
      call mem%dealloc(this%S, this%aux%n_ao, this%aux%n_ao)
      call mem%dealloc(this%Q, this%aux%n_ao, this%aux%n_ao)
      call mem%dealloc(this%Q_inverse, this%aux%n_ao, this%aux%n_ao)
!
      call this%aux%cleanup_shells()
!
   end subroutine cleanup_eri_ri
!
!
   subroutine construct_S(this, ao)
!!
!!    Construct S
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!!    Constructs the matrix
!!
!!       S_JK = (J|K)
!!
      implicit none
!
      class(eri_ri), intent(inout) :: this
      type(ao_tool), intent(in) :: ao
!
      integer :: J, K, w, x, w_red, x_red
!
      real(dp), dimension(this%aux%max_dim**2), target :: g_JK
      real(dp), dimension(:,:), pointer :: g_JK_p
!
      do J = 1, this%aux%n_shells
         do K = 1, this%aux%n_shells
!
            call ao%get_eri_2c(g_JK, this%aux%get_shell(J), this%aux%get_shell(K))
            g_JK_p(1 : this%aux%get_shell_length(J), 1 : this%aux%get_shell_length(K)) &
                        => g_JK(1 : this%aux%get_shell_length(J)*this%aux%get_shell_length(K))
!
            do w = this%aux%get_shell_first(J), this%aux%get_shell_last(J)
               do x = this%aux%get_shell_first(K), this%aux%get_shell_last(K)
!
                  w_red = w - this%aux%get_shell_first(J) + 1
                  x_red = x - this%aux%get_shell_first(K) + 1
!
                  this%S(w,x) = g_JK_p(w_red, x_red)
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine construct_S
!
!
   subroutine decompose_S(this)
!!
!!    Decompose S
!!    Written by Sarai D. Folketstad, Aug 2021
!!
!!    Construct Q by Cholesky decomposing S
!!
!!       S = QQ^T = PQ'(PQ')^T
!!
!!    Note that full_cholesky_decomposition is a pivoted routine,
!!    we must reorder the rows of Q to match the ordering
!!    of the integrals in libint
!
      use array_utilities, only: full_cholesky_decomposition
!
      implicit none
!
      class(eri_ri), intent(inout) :: this
!
      integer :: rank
      real(dp), dimension(:,:), allocatable :: Q_copy
      integer, dimension(:), allocatable :: pivots
!
      integer :: i, j
!
      call mem%alloc(Q_copy, this%aux%n_ao, this%aux%n_ao)
      call mem%alloc(pivots, this%aux%n_ao)
!
      call full_cholesky_decomposition(this%S, Q_copy, this%aux%n_ao, &
                                       rank, 1.0d-12, pivots)
!$omp parallel do private(i, j)
      do j = 1, this%aux%n_ao
         do i = 1, this%aux%n_ao
!
            this%Q(pivots(i), j) = Q_copy(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(Q_copy, this%aux%n_ao, this%aux%n_ao)
      call mem%dealloc(pivots, this%aux%n_ao)
!
   end subroutine decompose_S
!
!
   subroutine invert_Q(this)
!!
!!    Invert Q
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!
      use array_utilities, only: invert
!
      implicit none
!
      class(eri_ri), intent(inout) :: this
!
      call invert(this%Q_inverse, this%Q, this%aux%n_ao)
!
   end subroutine invert_Q
!
!
   subroutine print_summary(this)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      class(eri_ri) :: this
!
      call output%printf('n', '- Summary of RI approximation of &
                         &electronic repulsion integrals: ', fs='(/t3,a)')
      call output%printf('n', 'Auxilliary basis set: (a0)', &
                         chars=[this%aux%basis_set], fs='(/t6,a)')
      call output%printf('n', 'Dimension of auxilliary basis: (i0)', &
                         ints=[this%aux%n_ao], fs='(t6,a)')
!
   end subroutine print_summary
!
!
   pure function get_n_J_eri_ri(this) result(n_J)
!!
!!    Get n_J
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      class(eri_ri), intent(in) :: this
!
      integer :: n_J
!
      n_J = this%aux%n_ao
!
   end function get_n_J_eri_ri
!
!
end module eri_ri_class
