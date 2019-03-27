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
module mo_integral_tool_class
!
!!
!!    MO integral tool class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use file_class
   use disk_manager_class
   use memory_manager_class
   use timings_class
   use reordering
!
   implicit none
!
!  Class definition
!
   type :: mo_integral_tool
!
      logical, private :: cholesky_file      = .true.
      logical, private :: cholesky_t1_file   = .false.
!
      integer :: n_J
!
      type(file) :: cholesky_mo
      type(file) :: cholesky_mo_t1
!
      integer, private :: n_o
      integer, private :: n_v
      integer, private :: n_mo
!
      logical, private :: eri_t1_mem = .false.
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs
!
   contains
!
      procedure :: prepare                => prepare_mo_integral_tool
      procedure :: cleanup                => cleanup_mo_integral_tool
!
!     Read MO Cholesky vectors
!
      procedure :: read_cholesky          => read_cholesky_mo_integral_tool
!
!     Read/write/construct T1-transformed Cholesky vectors as well as T1-ERI construction
!
      procedure :: read_cholesky_t1       => read_cholesky_t1_mo_integral_tool
      procedure :: write_t1_cholesky      => write_t1_cholesky_mo_integral_tool
!
      procedure :: construct_g_pqrs_t1    => construct_g_pqrs_t1_mo_integral_tool
!
      procedure :: construct_cholesky_ij        => construct_cholesky_ij_mo_integral_tool
      procedure :: construct_cholesky_ab        => construct_cholesky_ab_mo_integral_tool
      procedure :: construct_cholesky_ai        => construct_cholesky_ai_mo_integral_tool
!
      procedure :: construct_cholesky_ij_c1     => construct_cholesky_ij_c1_mo_integral_tool
      procedure :: construct_cholesky_ab_c1     => construct_cholesky_ab_c1_mo_integral_tool
      procedure :: construct_cholesky_ai_a_c1   => construct_cholesky_ai_a_c1_mo_integral_tool
      procedure :: construct_cholesky_ai_i_c1   => construct_cholesky_ai_i_c1_mo_integral_tool
!
      procedure :: set_full_index         => set_full_index_mo_integral_tool
      procedure :: can_we_keep_g_pqrs_t1  => can_we_keep_g_pqrs_t1_mo_integral_tool
!
   end type mo_integral_tool
!
!
contains
!
!
   subroutine prepare_mo_integral_tool(integrals, n_J, n_o, n_v)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Initializes the integral tool. Note that the integral tool
!!    needs to know the number of occupied and virtual orbitals,
!!    which are also stored in the wavefunction.
!!
!!    n_o: number of occupied orbitals
!!    n_v: number of virtual orbitals
!!    n_J: number of Cholesky vectors
!!
      implicit none
!
      class(mo_integral_tool) :: integrals
!
      integer, intent(in) :: n_J
      integer, intent(in) :: n_o
      integer, intent(in) :: n_v
!
      integrals%n_J  = n_J
      integrals%n_o  = n_o
      integrals%n_v  = n_v
      integrals%n_mo = n_o + n_v
!
      call integrals%cholesky_mo%init('cholesky_mo_vectors', 'direct', 'unformatted', dp*n_J)
      call integrals%cholesky_mo_t1%init('cholesky_mo_t1_vectors', 'direct', 'unformatted', dp*n_J)
!
!     Initially MO cholesky on file, and not T1-transformed cholesky
!     nor full T1-ERI matrix
!
      integrals%cholesky_file      = .true.
      integrals%cholesky_t1_file   = .false.
      integrals%eri_t1_mem         = .false.
!
   end subroutine prepare_mo_integral_tool
!
!
   subroutine cleanup_mo_integral_tool(integrals)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(mo_integral_tool) :: integrals
!
      if (allocated(integrals%g_pqrs)) then
!
         call mem%dealloc(integrals%g_pqrs, integrals%n_mo, integrals%n_mo, &
                                          integrals%n_mo, integrals%n_mo)
!
      endif
!
   end subroutine cleanup_mo_integral_tool
!
!
   subroutine can_we_keep_g_pqrs_t1_mo_integral_tool(integrals)
!!
!!    Can we keep g_pqrs
!!    Written by Eirik F. Kjønstad, Jan 2019
!!
!!    This routine is called to check whether the T1-ERIs can be held in
!!    memory safely (< 20% of total available). If this is the case, the
!!    manager will keep a copy of g_pqrs in memory. When an integral is
!!    requested (e.g. g_abci), the integral will be copied from the g_pqrs
!!    copy instead of being constructed as g_abci = sum_J L_ab^J L_ci^J.
!!
!!    Note: the routine assumes that the T1-transformed Cholesky vectors
!!    have been placed on file.
!!
      implicit none
!
      class(mo_integral_tool) :: integrals
!
      real(dp), dimension(:,:,:), allocatable :: L_J_pq
!
      integer :: required_mem
!
      required_mem = (integrals%n_mo)**4 + (integrals%n_mo)**2*(integrals%n_J)
!
      if ((required_mem .lt. mem%get_available()/(5*dp)) .or. allocated(integrals%g_pqrs)) then
!
         integrals%eri_t1_mem = .true.
!
!        If not allocated, allocate copy of ERI matrix
!
         if (.not. allocated(integrals%g_pqrs)) then
!
            call mem%alloc(integrals%g_pqrs, integrals%n_mo, integrals%n_mo, &
                                          integrals%n_mo, integrals%n_mo)
!
         endif
!
!        Then construct it from the T1-Cholesky vectors
!
         call mem%alloc(L_J_pq, integrals%n_J, integrals%n_mo, integrals%n_mo)
!
         call integrals%read_cholesky_t1(L_J_pq, 1, integrals%n_mo, 1, integrals%n_mo)
!
         call dgemm('T','N',              &
                     integrals%n_mo**2,   &
                     integrals%n_mo**2,   &
                     integrals%n_J,       &
                     one,                 &
                     L_J_pq,              &
                     integrals%n_J,       &
                     L_J_pq,              &
                     integrals%n_J,       &
                     zero,                &
                     integrals%g_pqrs,    &
                     integrals%n_mo**2)
!
         call mem%dealloc(L_J_pq, integrals%n_J, integrals%n_mo, integrals%n_mo)
!
      else
!
         integrals%eri_t1_mem = .false.
!
      endif
!
   end subroutine can_we_keep_g_pqrs_t1_mo_integral_tool
!
!
   subroutine read_cholesky_mo_integral_tool(integrals, L_J_pq, first_p, last_p, first_q, last_q)
!!
!!    Read mo cholesky vectors
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Reads cholesky vectors L_pq_J for mo indices
!!    [first_p, first_p + dim_p - 1] and [first_q, first_q + dim_q - 1]
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      real(dp), dimension(integrals%n_J, last_p - first_p + 1, last_q - first_q + 1), intent(out) :: L_J_pq
!
      integer :: p, q, pq_rec
!
      call disk%open_file(integrals%cholesky_mo, 'read')
!
      do q = 1, last_q - first_q + 1
         do p = 1, last_p - first_p + 1
!
            pq_rec = (max(p + first_p - 1, q + first_q - 1)*(max(p + first_p - 1, q + first_q - 1)-3)/2) &
                         + (p + first_p - 1) + (q + first_q - 1)
!
            read(integrals%cholesky_mo%unit, rec=pq_rec) L_J_pq(:, p, q)
!
         enddo
      enddo
!
      call disk%close_file(integrals%cholesky_mo)
!
   end subroutine read_cholesky_mo_integral_tool
!
!
   subroutine read_cholesky_t1_mo_integral_tool(integrals, L_J_pq, first_p, last_p, first_q, last_q)
!!
!!    Read mo t1 cholesky vectors
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Reads T1-transformed cholesky vectors L_pq_J for mo indices
!!    [first_p, last_p] and [first_q, last_q]
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      real(dp), dimension(integrals%n_J, last_p - first_p + 1, last_q - first_q + 1), intent(out) :: L_J_pq
!
      integer :: p, q, pq_rec
!
      call disk%open_file(integrals%cholesky_mo_t1, 'read')
!
      do q = 1, last_q - first_q + 1
         do p = 1, last_p - first_p + 1
!
            pq_rec = integrals%n_mo*(q + first_q - 2) + p + first_p - 1
!
            read(integrals%cholesky_mo_t1%unit, rec=pq_rec) L_J_pq(:, p, q)
!
         enddo
      enddo
!
      call disk%close_file(integrals%cholesky_mo_t1)
!
   end subroutine read_cholesky_t1_mo_integral_tool
!
!
   subroutine construct_g_pqrs_t1_mo_integral_tool(integrals, g_pqrs, first_p, last_p, first_q, last_q, &
                                                                     first_r, last_r, first_s, last_s)
!!
!!    Construct g_pqrs T1
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
!!    Assumes that the integrals are in memory or that the Cholesky vectors
!!    are stored in the T1 transformed basis.
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
      integer, intent(in) :: first_r, last_r
      integer, intent(in) :: first_s, last_s
!
      real(dp), dimension(last_p-first_p+1,last_q-first_q+1,last_r-first_r+1,last_s-first_s+1), intent(out) :: g_pqrs
!
      integer :: p, q, r, s
      integer :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:,:), allocatable :: L_J_pq, L_J_rs
!
      if (.not. integrals%cholesky_t1_file) call output%error_msg('tried to construct T1-tranformed g_pqrs, but ' &
                                                         // 'the T1-transformed Cholesky vectors are not on file!')
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
      dim_r = last_r - first_r + 1
      dim_s = last_s - first_s + 1
!
      if (integrals%eri_t1_mem) then
!
!$omp parallel do private(s, r, q, p)
         do s = 1, dim_s
            do r = 1, dim_r
               do q = 1, dim_q
                  do p = 1, dim_p
!
                     g_pqrs(p,q,r,s) = integrals%g_pqrs(p + first_p - 1,   &
                                                         q + first_q - 1,  &
                                                         r + first_r - 1,  &
                                                         s + first_s - 1)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
      else
!
         call mem%alloc(L_J_pq, integrals%n_J, dim_p, dim_q)
         call mem%alloc(L_J_rs, integrals%n_J, dim_r, dim_s)
!
         call integrals%read_cholesky_t1(L_J_pq, first_p, last_p, first_q, last_q)
         call integrals%read_cholesky_t1(L_J_rs, first_r, last_r, first_s, last_s)
!
         call dgemm('T', 'N', &
                     dim_p*dim_q, &
                     dim_r*dim_s, &
                     integrals%n_J, &
                     one, &
                     L_J_pq, &
                     integrals%n_J, &
                     L_J_rs, &
                     integrals%n_J, &
                     zero, &
                     g_pqrs, &
                     dim_p*dim_q)
!
         call mem%dealloc(L_J_pq, integrals%n_J, dim_p, dim_q)
         call mem%dealloc(L_J_rs, integrals%n_J, dim_r, dim_s)
!
      endif
!
   end subroutine construct_g_pqrs_t1_mo_integral_tool
!
!
   subroutine construct_cholesky_ij_mo_integral_tool(integrals, L_ij_J, t1, first_i, last_i, first_j, last_j)
!!
!!    Construct Cholesky ij
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Computes
!!
!!       L_ij_J_T1 = L_ij_J + sum_a t_aj L_ia_J,
!!
!!    and saves the result in L_ij_J.
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_j, last_j
!
      real(dp), dimension(last_i - first_i + 1, last_j - first_j + 1, integrals%n_J), intent(inout) :: L_ij_J
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: t1
!
      real(dp), dimension(:,:,:), allocatable :: L_iJ_j_term
      real(dp), dimension(:,:,:), allocatable :: L_iJ_a
      real(dp), dimension(:,:,:), allocatable :: L_J_ia, L_J_ij
!
      integer :: full_first_i, full_last_i
      integer :: full_first_j, full_last_j
!
      integer :: i_length, j_length
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_j, 'f', 'o', first_j)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_j, 'l', 'o', last_j)
!
      i_length = full_last_i - full_first_i + 1
      j_length = full_last_j - full_first_j + 1
!
!     Read the untransformed Cholesky vectors
!
      call mem%alloc(L_J_ij, integrals%n_J, i_length, j_length)
      call integrals%read_cholesky(L_J_ij, full_first_i, full_last_i, full_first_j, full_last_j)
      call sort_12_to_21(L_J_ij, L_ij_J, integrals%n_J, i_length*j_length)
      call mem%dealloc(L_J_ij, integrals%n_J, i_length, j_length)
!
      call mem%alloc(L_J_ia, integrals%n_J, i_length, integrals%n_v)
      call integrals%read_cholesky(L_J_ia, full_first_i, full_last_i, integrals%n_o + 1, integrals%n_mo)
!
!     Compute and add t1-transformed term, L_iJ_j sum_a t_aj L_ia_J
!
      call mem%alloc(L_iJ_a, i_length, integrals%n_J, integrals%n_v)
      call sort_123_to_213(L_J_ia, L_iJ_a, integrals%n_J, i_length, integrals%n_v)
      call mem%dealloc(L_J_ia, integrals%n_J, i_length, integrals%n_v)
!
      call mem%alloc(L_iJ_j_term, i_length, integrals%n_J, j_length)
!
      call dgemm('N','N',                   &
                  i_length*(integrals%n_J), &
                  j_length,                 &
                  integrals%n_v,            &
                  one,                      &
                  L_iJ_a,                   &
                  i_length*(integrals%n_J), &
                  t1(1, first_j),           &
                  integrals%n_v,            &
                  zero,                     &
                  L_iJ_j_term,              &
                  i_length*(integrals%n_J))
!
     call add_132_to_123(one, L_iJ_j_term, L_ij_J, i_length, j_length, integrals%n_J)
!
     call mem%dealloc(L_iJ_j_term, i_length, integrals%n_J, j_length)
     call mem%dealloc(L_iJ_a, i_length, integrals%n_J, integrals%n_v)
!
   end subroutine construct_cholesky_ij_mo_integral_tool
!
!
   subroutine construct_cholesky_ab_mo_integral_tool(integrals, L_ab_J, t1, first_a, last_a, first_b, last_b)
!!
!!    Construct Cholesky ab
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Computes
!!
!!       L_ab_J_T1 = L_ab_J - sum_i t_ai L_ib_J
!!
!!    and saves the result in L_ab_J. Note that batching is handled
!!    outside this routine, not inside.
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: t1
!
      integer, intent(in) :: first_a, last_a
      integer, intent(in) :: first_b, last_b
!
      real(dp), dimension(last_a - first_a + 1, last_b - first_b + 1, integrals%n_J) :: L_ab_J
!
      integer :: full_first_a, full_last_a
      integer :: full_first_b, full_last_b
!
      real(dp), dimension(:,:,:), allocatable :: L_ib_J, L_J_ib, L_J_ab
!
      integer :: b_length, a_length
!
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
      call integrals%set_full_index(full_first_b, 'f', 'v', first_b)
!
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
      call integrals%set_full_index(full_last_b, 'l', 'v', last_b)
!
      a_length = full_last_a - full_first_a + 1
      b_length = full_last_b - full_first_b + 1
!
!     Set L_ib_J = L_ib^J and L_ab_J^T1 = L_ab_J
!
      call mem%alloc(L_J_ib, integrals%n_J, integrals%n_o, b_length)
      call integrals%read_cholesky(L_J_ib, 1, integrals%n_o, full_first_b, full_last_b)
!
      call mem%alloc(L_ib_J, integrals%n_o, b_length, integrals%n_J)
      call sort_12_to_21(L_J_ib, L_ib_J, integrals%n_J, integrals%n_o*b_length)
      call mem%dealloc(L_J_ib, integrals%n_J, integrals%n_o, b_length)
!
      call mem%alloc(L_J_ab, integrals%n_J, a_length, b_length)
      call integrals%read_cholesky(L_J_ab, full_first_a, full_last_a, full_first_b, full_last_b)
      call sort_12_to_21(L_J_ab, L_ab_J, integrals%n_J, a_length*b_length)
      call mem%dealloc(L_J_ab, integrals%n_J, a_length, b_length)
!
!     Calculate and add t1-transformed term, - sum_i t_ai L_ib_J
!
      call dgemm('N','N',                   &
                  a_length,                 &
                  b_length*(integrals%n_J), &
                  integrals%n_o,            &
                  -one,                     &
                  t1(first_a, 1),           & ! t_a_i
                  integrals%n_v,            &
                  L_ib_J,                   & ! L_i_bJ
                  integrals%n_o,            &
                  one,                      &
                  L_ab_J,                   & ! L_a_bj
                  a_length)
!
      call mem%dealloc(L_ib_J, integrals%n_o, b_length, integrals%n_J)
!
   end subroutine construct_cholesky_ab_mo_integral_tool
!
!
   subroutine construct_cholesky_ai_mo_integral_tool(integrals, L_ai_J, t1, first_a, last_a, first_i, last_i)
!!
!!    Construct Cholesky ai
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Computes the T1-transformed cholesky vector
!!
!!       L_ai_J_T1 = L_ai_J - sum_j t_aj*L_ji_J
!!                          + sum_b t_bi*L_ab_J
!!                          - sum_bj t_aj*t_bi*L_jb_J
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_a, last_a
!
      real(dp), dimension(last_a - first_a + 1, last_i - first_i + 1, integrals%n_J) :: L_ai_J
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: t1
!
      real(dp), dimension(:, :, :), allocatable :: L_J_ji, L_ji_J, X_j_iJ, L_J_ab, X_J_ai, L_J_ai, X_Jj_i, L_J_jb
!
      integer :: full_first_a, full_last_a, length_a
      integer :: full_first_i, full_last_i, length_i
!
      type(batching_index) :: batch_j
!
      integer :: req0, req1, current_j_batch
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
      length_i = full_last_i - full_first_i + 1
      length_a = full_last_a - full_first_a + 1
!
!     :: Term 1: L_ai_J - sum_bj t_aj*t_bi*L_jb_J
!
      call mem%alloc(L_J_ai, integrals%n_J, length_a, length_i)
      call integrals%read_cholesky(L_J_ai, full_first_a, full_last_a, full_first_i, full_last_i)
      call sort_12_to_21(L_J_ai, L_ai_J, integrals%n_J, length_a*length_i)
      call mem%dealloc(L_J_ai, integrals%n_J, length_a, length_i)
!
      call batch_j%init(integrals%n_o)
!
      req0 = 0
      req1 = (integrals%n_o)*(integrals%n_J) & ! X_i_jJ
            + (integrals%n_v)*(integrals%n_J)  ! L_bj_J
!
      call mem%batch_setup(batch_j, req0, req1)
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_J_jb, integrals%n_J, batch_j%length, integrals%n_v)
         call integrals%read_cholesky(L_J_jb, batch_j%first, batch_j%last, (integrals%n_o) + 1, (integrals%n_mo))
!
         call mem%alloc(X_Jj_i, integrals%n_J, batch_j%length, length_i)
!
         call dgemm('N', 'N',                         &
                  (batch_j%length)*(integrals%n_J),   &
                  length_i,                           &
                  integrals%n_v,                      &
                  one,                                &
                  L_J_jb,                             &   ! L_Jj_b
                  (batch_j%length)*(integrals%n_J),   &
                  t1(1, first_i),                     &   ! t_b_i
                  integrals%n_v,                      &
                  zero,                               &
                  X_Jj_i,                             &
                  (batch_j%length)*(integrals%n_J))
!
         call mem%dealloc(L_J_jb, integrals%n_J, batch_j%length, integrals%n_v)
!
         call mem%alloc(X_j_iJ, batch_j%length, length_i, integrals%n_J)
!
         call sort_12_to_21(X_Jj_i, X_j_iJ, (integrals%n_J), batch_j%length*length_i)
!
         call mem%dealloc(X_Jj_i, integrals%n_J, batch_j%length, length_i)
!
         call dgemm('N', 'N',                   &
                  length_a,                     &
                  (length_i)*(integrals%n_J),   &
                  batch_j%length,               &
                  -one,                         &
                  t1(first_a, batch_j%first),   & ! t_a_j
                  integrals%n_v,                &
                  X_j_iJ,                       &
                  batch_j%length,               &
                  one,                          &
                  L_ai_J,                       & ! L_a_iJ
                  length_a)
!
         call mem%dealloc(X_j_iJ, batch_j%length, length_i, integrals%n_J)
!
      enddo
!
!     :: Term 2: L_ai_J - sum_j t_aj*L_ji_J
!
      call mem%alloc(L_J_ji, integrals%n_J, integrals%n_o, length_i)
      call integrals%read_cholesky(L_J_ji, 1, (integrals%n_o), full_first_i, full_last_i)
!
      call mem%alloc(L_ji_J, integrals%n_o, length_i, integrals%n_J)
      call sort_12_to_21(L_J_ji, L_ji_J, integrals%n_J, integrals%n_o*length_i)
      call mem%dealloc(L_J_ji, integrals%n_J, integrals%n_o, length_i)
!
       call dgemm('N', 'N',                     &
                   length_a,                    &
                   (length_i)*(integrals%n_J),  &
                   integrals%n_o,               &
                   -one,                        &
                   t1(first_a, 1),              & ! t_a_j
                   integrals%n_v,               &
                   L_ji_J,                      & ! L_j_iJ
                   integrals%n_o,               &
                   one,                         &
                   L_ai_J,                      & ! L_a_iJ
                   length_a)
!
      call mem%dealloc(L_ji_J, integrals%n_o, length_i, integrals%n_J)
!
!     :: Term 3: L_ai_J + sum_b t_bi*L_ab_J
!
      call mem%alloc(L_J_ab, integrals%n_J, length_a, integrals%n_v)
!
      call integrals%read_cholesky(L_J_ab, &
                                       first_a + integrals%n_o, last_a + (integrals%n_o), &
                                       integrals%n_o + 1, integrals%n_mo)
!
!     X_Ja,i = L_Ja,b t_b,i
!
      call mem%alloc(X_J_ai, integrals%n_J, length_a, length_i)
!
      call dgemm('N','N', &
                  integrals%n_J*length_a, &
                  length_i, &
                  integrals%n_v, &
                  one, &
                  L_J_ab, &
                  integrals%n_J*length_a, &
                  t1(1,first_i), &
                  integrals%n_v, &
                  zero, &
                  X_J_ai, &
                  integrals%n_J*length_a)
!
      call mem%dealloc(L_J_ab, integrals%n_J, length_a, integrals%n_v)
!
      call add_21_to_12(one, X_J_ai, L_ai_J, length_a*length_i, integrals%n_J)
!
      call mem%dealloc(X_J_ai, integrals%n_J, length_a, length_i)
!
   end subroutine construct_cholesky_ai_mo_integral_tool
!
!
   subroutine construct_cholesky_ij_c1_mo_integral_tool(integrals, L_J_ij_c1, c_aj, first_i, last_i, first_j, last_j)
!!
!!    Construct the C1-transformed Cholesky Vector ij from the T1-transformed Cholesky Vector L_J_ia
!!
!!    Based on construct_cholesky_ij_mo_integral_tool written by Sarai D. Folkestad and Eirik F. Kjønstad
!!    Modified for the c1-transformation by Alexander Paul, Feb 2019
!!
!!    j is the transformed index
!!
!!    Note: the routine assumes that the T1-transformed Cholesky vectors have been placed on file.
!!
!!    Computes
!!
!!       L_ij_J_c1= sum_a L_ia_J_T1 c_aj ,
!!
!!    and saves the result in L_J_ij_c1 (ordering to be consistent with read_cholesky_t1)
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_j, last_j
!
      real(dp), dimension(integrals%n_J, last_i - first_i + 1, last_j - first_j + 1), intent(inout) :: L_J_ij_c1
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: c_aj
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ia
!
      integer :: full_first_i, full_last_i
      integer :: full_first_j, full_last_j
!
      integer :: i_length, j_length
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_j, 'f', 'o', first_j)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_j, 'l', 'o', last_j)
!
      i_length = full_last_i - full_first_i + 1
      j_length = full_last_j - full_first_j + 1
!
!     Read the t1-transformed Cholesky vectors
!
      call mem%alloc(L_J_ia, integrals%n_J, i_length, integrals%n_v)
      call integrals%read_cholesky_t1(L_J_ia, full_first_i, full_last_i, integrals%n_o + 1, integrals%n_mo)
!
!     Compute and c1-transformed term, L_J_ij = sum_a c_aj L_J_ia
!
      call dgemm('N','N',                    &
                  i_length*(integrals%n_J),  &
                  j_length,                  &
                  integrals%n_v,             &
                  one,                       &
                  L_J_ia,                    & ! L_Ji_a
                  i_length*(integrals%n_J),  &
                  c_aj(1, first_j),          & ! c_a_j
                  integrals%n_v,             &
                  zero,                      &
                  L_J_ij_c1,                 & ! L_Ji_j
                  i_length*(integrals%n_J))
!
      call mem%dealloc(L_J_ia, integrals%n_J, i_length, integrals%n_v)
!
   end subroutine construct_cholesky_ij_c1_mo_integral_tool
!
!
   subroutine construct_cholesky_ab_c1_mo_integral_tool(integrals, L_J_ab_c1, c_ai, first_a, last_a, first_b, last_b)
!!
!!    Construct the C1-transformed Cholesky Vector ab from the T1-transformed Cholesky Vector
!!
!!    Based on construct_cholesky_ab_mo_integral_tool written by Sarai D. Folkestad and Eirik F. Kjønstad
!!    Modified for the c1-transformation by Alexander Paul, Feb 2019
!!
!!    a is the transformed index
!!
!!    Note: the routine assumes that the T1-transformed Cholesky vectors have been placed on file.
!!
!!    Computes
!!
!!       L_ab_J_c1= - sum_i c_ai L_ib_J_T1 ,
!!
!!    and returns the result in L_J_ab_c1 (ordering to be consistent with read_cholesky_t1)
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: c_ai
!
      integer, intent(in) :: first_a, last_a
      integer, intent(in) :: first_b, last_b
!
      real(dp), dimension(integrals%n_J, last_a - first_a + 1, last_b - first_b + 1) :: L_J_ab_c1
!
      integer :: full_first_a, full_last_a
      integer :: full_first_b, full_last_b
!
      real(dp), dimension(:,:,:), allocatable :: L_ib_J, L_J_ib, L_ab_J_c1
!
      integer :: b_length, a_length
!
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
      call integrals%set_full_index(full_first_b, 'f', 'v', first_b)
!
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
      call integrals%set_full_index(full_last_b, 'l', 'v', last_b)
!
      a_length = full_last_a - full_first_a + 1
      b_length = full_last_b - full_first_b + 1
!
!     Read t1-transformed L_J_ib and resort
!
      call mem%alloc(L_J_ib, integrals%n_J, integrals%n_o, b_length)
      call integrals%read_cholesky_t1(L_J_ib, 1, integrals%n_o, full_first_b, full_last_b)
!
      call mem%alloc(L_ib_J, integrals%n_o, b_length, integrals%n_J)
      call sort_12_to_21(L_J_ib, L_ib_J, integrals%n_J, integrals%n_o*b_length)
      call mem%dealloc(L_J_ib, integrals%n_J, integrals%n_o, b_length)
!
!     Calculate and add c1-transformed L_ab_J_c1 = - sum_i c_ai L_ib_J
!
      call mem%alloc(L_ab_J_c1, a_length, b_length, integrals%n_J)
!
      call dgemm('N','N',                   &
                  a_length,                 &
                  b_length*(integrals%n_J), &
                  integrals%n_o,            &
                  -one,                     &
                  c_ai(first_a, 1),         & ! c_a_i
                  integrals%n_v,            &
                  L_ib_J,                   & ! L_i_bJ
                  integrals%n_o,            &
                  zero,                     &
                  L_ab_J_c1,                & ! L_a_bj
                  a_length)
!
      call mem%dealloc(L_ib_J, integrals%n_o, b_length, integrals%n_J)
!
      call sort_12_to_21(L_ab_J_c1, L_J_ab_c1, a_length*b_length, integrals%n_J)
!
      call mem%dealloc(L_ab_J_c1, a_length, b_length, integrals%n_J)
!
   end subroutine construct_cholesky_ab_c1_mo_integral_tool
!
!
   subroutine construct_cholesky_ai_a_c1_mo_integral_tool(integrals, L_J_ai_c1, c_aj, first_a, last_a, first_i, last_i)
!!
!!    Construct the C1-transformed Cholesky Vector ai from the T1-transformed Cholesky Vector
!!
!!    Based on construct_cholesky_ai_mo_integral_tool (Term 3) written by Sarai D. Folkestad and Eirik F. Kjønstad
!!    Modified for the c1-transformation by Alexander Paul, Feb 2019
!!
!!    a is the transformed index
!!
!!    Note: the routine assumes that the T1-transformed Cholesky vectors
!!    have been placed on file.
!!
!!    Computes
!!
!!       L_ai_J_c1= - sum_j c_aj L_ji_J_T1 ,
!!
!!    and returns the result in L_J_ai_c1 (ordering to be consistent with read_cholesky_t1)
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_a, last_a
!
      real(dp), dimension(last_a - first_a + 1, last_i - first_i + 1, integrals%n_J) :: L_J_ai_c1
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: c_aj
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ji, L_ji_J, L_ai_J_c1
!
      integer :: full_first_a, full_last_a, length_a
      integer :: full_first_i, full_last_i, length_i
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
      length_i = full_last_i - full_first_i + 1
      length_a = full_last_a - full_first_a + 1
!
!     Read t1-transformed L_J_ji and resort
!
      call mem%alloc(L_J_ji, integrals%n_J, integrals%n_o, length_i)
      call integrals%read_cholesky_t1(L_J_ji, 1, integrals%n_o, full_first_i, full_last_i)
!
      call mem%alloc(L_ji_J, integrals%n_o, length_i, integrals%n_J)
      call sort_12_to_21(L_J_ji, L_ji_J, integrals%n_J, integrals%n_o*length_i)
      call mem%dealloc(L_J_ji, integrals%n_J, integrals%n_o, length_i)
!
      call mem%alloc(L_ai_J_c1, length_a, length_i, integrals%n_J)
!
      call dgemm('N', 'N',                     &
                  length_a,                    &
                  (length_i)*(integrals%n_J),  &
                  integrals%n_o,               &
                  -one,                        &
                  c_aj(first_a, 1),            & ! c_a_j
                  integrals%n_v,               &
                  L_ji_J,                      & ! L_j_iJ
                  integrals%n_o,               &
                  zero,                        &
                  L_ai_J_c1,                   & ! L_a_iJ
                  length_a)
!
      call mem%dealloc(L_ji_J, integrals%n_o, length_i, integrals%n_J)
!
      call sort_12_to_21(L_ai_J_c1, L_J_ai_c1, length_a*length_i, integrals%n_J)
!
      call mem%dealloc(L_ai_J_c1, length_a, length_i, integrals%n_J)
!
   end subroutine construct_cholesky_ai_a_c1_mo_integral_tool
!
!
   subroutine construct_cholesky_ai_i_c1_mo_integral_tool(integrals, L_J_ai_c1, c_bi, first_a, last_a, first_i, last_i)
!!
!!    Construct the C1-transformed Cholesky Vector ai from the T1-transformed Cholesky Vector
!!
!!    Based on construct_cholesky_ai_mo_integral_tool (Term 2) written by Sarai D. Folkestad and Eirik F. Kjønstad
!!    Modified for the c1-transformation by Alexander Paul, Feb 2019
!!
!!    i is the transformed index
!!
!!    Note: the routine assumes that the T1-transformed Cholesky vectors
!!    have been placed on file.
!!
!!    Computes
!!
!!       L_ai'_J_c1= sum_b L_ab_J_T1 c_bi,
!!
!!    and returns the result in L_J_ai_c1 (ordering to be consistent with read_cholesky_t1)
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(in) :: first_i, last_i
      integer, intent(in) :: first_a, last_a
!
      real(dp), dimension(last_a - first_a + 1, last_i - first_i + 1, integrals%n_J) :: L_J_ai_c1
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: c_bi
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ab
!
      integer :: full_first_a, full_last_a, length_a
      integer :: full_first_i, full_last_i, length_i
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
      length_i = full_last_i - full_first_i + 1
      length_a = full_last_a - full_first_a + 1
!
!     Read t1-transformed L_J_ab
!
      call mem%alloc(L_J_ab, integrals%n_J, length_a, integrals%n_v)
      call integrals%read_cholesky_t1(L_J_ab, &
                                       first_a + integrals%n_o, last_a + (integrals%n_o), &
                                       integrals%n_o + 1, integrals%n_mo)
!
       call dgemm('N', 'N',               &
                  integrals%n_J*length_a, &
                  length_i,               &
                  integrals%n_v,          &
                  one,                    &
                  L_J_ab,                 & ! L_Ja_b
                  integrals%n_J*length_a, &
                  c_bi(1, first_i),       & ! c_b_i
                  integrals%n_v,          &
                  zero,                   &
                  L_J_ai_c1,              & ! L_Ja_i
                  integrals%n_J*length_a)
!
      call mem%dealloc(L_J_ab, integrals%n_J, length_a, integrals%n_v)
!
end subroutine construct_cholesky_ai_i_c1_mo_integral_tool
!
!
   subroutine set_full_index_mo_integral_tool(integrals, ind, pos, orb_space, red_ind)
!!
!!    Set full index
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Sets the full space index "ind" based on the provided orbital space ('o' or 'v'
!!    for occupied or virtual). Prepares the index "ind" to call the routine to
!!    read cholesky vector. The "pos" specifies whether it is the first or last orbital
!!    index ('f' or 'l').
!!
      implicit none
!
      class(mo_integral_tool), intent(in) :: integrals
!
      integer, intent(inout) :: ind
!
      character(len=1), intent(in) :: orb_space
      character(len=1), intent(in) :: pos
!
      integer, optional, intent(in) :: red_ind
!
      if (present(red_ind)) then
!
         if (trim(orb_space) == 'o') then
!
            ind = red_ind
!
         elseif (trim(orb_space) == 'v') then
!
            ind = integrals%n_o + red_ind
!
         else
!
            call output%error_msg('did not recognize orbital space' // orb_space // ' in integral tool')
!
         endif
!
      else
!
         if (trim(orb_space) == 'o') then
!
            if (trim(pos) == 'f') then ! First index
!
               ind = 1
!
            elseif (trim(pos) == 'l') then ! Last index
!
               ind = integrals%n_o
!
            else
!
               call output%error_msg('did not recognize position' // pos // ' in integral tool')
!
            endif
!
         elseif (trim(orb_space) == 'v') then
!
            if (trim(pos) == 'f') then ! First index
!
               ind = integrals%n_o + 1
!
            elseif (trim(pos) == 'l') then ! Last index
!
               ind = integrals%n_o + integrals%n_v
!
            else
!
               call output%error_msg('did not recognize position' // pos // ' in integral tool')
!
            endif
!
         else
!
            call output%error_msg('did not recognize orbital space' // orb_space // ' in integral tool')
!
         endif
!
      endif
!
   end subroutine set_full_index_mo_integral_tool
!
!
   subroutine write_t1_cholesky_mo_integral_tool(integrals, t1)
!!
!!    Write T1-transformed Cholesky vectors to file
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Eirik F. Kjønstad, Mar 2019: Modifications to write consequtively to file.
!!
      implicit none
!
      class(mo_integral_tool), intent(inout) :: integrals
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) ::t1
!
      real(dp), dimension(:,:,:), allocatable :: L_ij_J, L_ai_J, L_ab_J
      real(dp), dimension(:,:,:), allocatable :: L_J_ia, L_J_ab, L_J_ij, L_J_ai
!
      integer :: ij_rec, i, j, ai_rec, ia_rec, ab_rec, a, b
!
      integer :: req0, req1, req1_a, req1_i, req2, current_i_batch, current_a_batch, current_b_batch
!
      type(batching_index) :: batch_i, batch_a, batch_b
!
      type(timings) :: write_t1_cholesky_timer
!
      call write_t1_cholesky_timer%init('transform and write t1 cholesky to file')
      call write_t1_cholesky_timer%start()
!
      call disk%open_file(integrals%cholesky_mo_t1, 'write')
!
!     occupied-occupied block
!
      call batch_i%init(integrals%n_o)
!
      req0 = 0
!
      req1 = 2*(integrals%n_o)*(integrals%n_J)    &   ! 2 x L_ij^J
            + 2*(integrals%n_v)*(integrals%n_J)       ! 2 x L_ia^J
!
      call mem%batch_setup(batch_i, req0, req1)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call mem%alloc(L_ij_J, batch_i%length, integrals%n_o, integrals%n_J)
         call integrals%construct_cholesky_ij(L_ij_J, t1, batch_i%first, batch_i%last, 1, integrals%n_o)
!
         call mem%alloc(L_J_ij, integrals%n_J, batch_i%length, integrals%n_o)
         call sort_12_to_21(L_ij_J, L_J_ij, integrals%n_o*batch_i%length, integrals%n_J)
         call mem%dealloc(L_ij_J, batch_i%length, integrals%n_o, integrals%n_J)
!
         do j = 1, integrals%n_o
            do i = 1, batch_i%length
!
               ij_rec = integrals%n_mo*(j - 1) + i + batch_i%first - 1
!
               write(integrals%cholesky_mo_t1%unit, rec=ij_rec) L_J_ij(:, i, j)
!
            enddo
         enddo
!
         call mem%dealloc(L_J_ij, integrals%n_J, batch_i%length, integrals%n_o)
!
      enddo
!
!     occupied-virtual block
!
      call batch_a%init(integrals%n_v)
!
      req0 = 0
!
      req1 = (integrals%n_o)*(integrals%n_J) ! L_ia^J
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(L_J_ia, integrals%n_J, integrals%n_o, batch_a%length)
!
         call integrals%read_cholesky(L_J_ia, 1, integrals%n_o, integrals%n_o + batch_a%first, integrals%n_o + batch_a%last)
!
         do a = 1, batch_a%length
            do i = 1, integrals%n_o
!
               ia_rec = integrals%n_mo*(a + batch_a%first + integrals%n_o - 2) + i
!
               write(integrals%cholesky_mo_t1%unit, rec=ia_rec) L_J_ia(:, i, a)
!
            enddo
         enddo
!
         call mem%dealloc(L_J_ia, integrals%n_J, integrals%n_o, batch_a%length)
!
      enddo
!
!     virtual-occupied block
!
      req0 = (integrals%n_v)*(integrals%n_J)   ! Uncontrollable from outside (L_jb^J, batching over j inside)
!
      req1_a = (integrals%n_v)*(integrals%n_J)   ! L_ab^J
      req1_i = 2*(integrals%n_o)*(integrals%n_J) ! L_ji^J
!
      req2 = 2*(integrals%n_J) ! 2 x L_ai^J
!
      call batch_i%init(integrals%n_o)
      call batch_a%init(integrals%n_v)
!
      call mem%batch_setup(batch_a, batch_i, req0, req1_a, req1_i, req2)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_i_batch = 1, batch_i%num_batches
!
            call batch_i%determine_limits(current_i_batch)
!
            call mem%alloc(L_ai_J, batch_a%length, batch_i%length, integrals%n_J)
!
            call integrals%construct_cholesky_ai(L_ai_J, t1, batch_a%first, batch_a%last, batch_i%first, batch_i%last)
!
            call mem%alloc(L_J_ai, integrals%n_J, batch_a%length, batch_i%length)
            call sort_12_to_21(L_ai_J, L_J_ai, batch_a%length*batch_i%length, integrals%n_J)
            call mem%dealloc(L_ai_J, batch_a%length, batch_i%length, integrals%n_J)
!
            do i = 1, batch_i%length
               do a = 1, batch_a%length
!
                  ai_rec = integrals%n_mo*(i + batch_i%first - 2) + a + batch_a%first - 1 + integrals%n_o
!
                  write(integrals%cholesky_mo_t1%unit, rec=ai_rec) L_J_ai(:, a, i)
!
               enddo
            enddo
!
            call mem%dealloc(L_J_ai, integrals%n_J, batch_a%length, batch_i%length)
!
         enddo
      enddo
!
!     virtual-virtual block
!
      call batch_b%init(integrals%n_v)
!
      req0 = 0
!
      req1 = 2*(integrals%n_v)*(integrals%n_J) &   ! 2 x L_ab^J
            + 2*(integrals%n_o)*(integrals%n_J)      ! L_ib^J
!
      call mem%batch_setup(batch_b, req0, req1)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         call mem%alloc(L_ab_J, integrals%n_v, batch_b%length, integrals%n_J)
!
         call integrals%construct_cholesky_ab(L_ab_J, t1, 1, integrals%n_v, batch_b%first, batch_b%last)
!
         call mem%alloc(L_J_ab, integrals%n_J, integrals%n_v, batch_b%length)
         call sort_12_to_21(L_ab_J, L_J_ab, integrals%n_v*batch_b%length, integrals%n_J)
         call mem%dealloc(L_ab_J, integrals%n_v, batch_b%length, integrals%n_J)
!
         do b = 1, batch_b%length
            do a = 1, integrals%n_v
!
               ab_rec = integrals%n_mo*((integrals%n_o + b + batch_b%first) - 2) + (a + integrals%n_o)
!
               write(integrals%cholesky_mo_t1%unit, rec=ab_rec) L_J_ab(:, a, b)
!
            enddo
         enddo
!
         call mem%dealloc(L_J_ab, integrals%n_J, integrals%n_v, batch_b%length)
!
      enddo
!
      integrals%cholesky_t1_file = .true.
!
      call disk%close_file(integrals%cholesky_mo_t1)
!
      call write_t1_cholesky_timer%freeze()
      call write_t1_cholesky_timer%switch_off()
!
   end subroutine write_t1_cholesky_mo_integral_tool
!
!
end module mo_integral_tool_class
