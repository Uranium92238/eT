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
module cholesky_orbital_tool_class
!
!!
!!    Cholesky orbital tool class
!!    Written by Sarai D. Folkestad, Mar 2021
!!
!!    Handles decomposition of the density to generate Cholesky orbitals
!!
!
   use parameters
   use memory_manager_class, only: mem
   use global_out,           only: output
!
   implicit none
!
   type :: cholesky_orbital_tool
!
      real(dp) :: threshold
!
      real(dp), dimension(:,:), allocatable :: D
!
      integer :: n_ao
!
      integer :: initial_rank, rank
!
      contains
!
         procedure :: restricted_decomposition &
                   => restricted_decomposition_cholesky_orbital_tool
!
         procedure :: full_decomposition &
                   => full_decomposition_cholesky_orbital_tool
!
         procedure :: set_density_from_density &
                   => set_density_from_density_cholesky_orbital_tool
!
         procedure :: set_density_from_orbitals &
                   => set_density_from_orbitals_cholesky_orbital_tool
!
         procedure :: initialize_density &
                   => initialize_density_cholesky_orbital_tool
!
         procedure :: destruct_density &
                   => destruct_density_cholesky_orbital_tool
!
         procedure :: cleanup &
                   => cleanup_cholesky_orbital_tool
!
         procedure, nopass, private :: cholesky_decomposition_limited_diagonal
         procedure, nopass, private :: standard_cholesky_decomposition
!
   end type  cholesky_orbital_tool
!
   interface  cholesky_orbital_tool
!
      procedure :: new_cholesky_orbital_tool
!
   end interface  cholesky_orbital_tool
!
contains
!
!
   pure function new_cholesky_orbital_tool(n_ao, threshold) result(this)
!!
!!    New cholesky orbital tool
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      type(cholesky_orbital_tool) :: this
!
      integer, intent(in)  :: n_ao
      real(dp), intent(in) :: threshold
!
      this%n_ao = n_ao
      this%threshold = threshold
!
   end function new_cholesky_orbital_tool
!
!
   subroutine initialize_density_cholesky_orbital_tool(this)
!!
!!    Initialize density
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(cholesky_orbital_tool), intent(inout)  :: this
!
      call mem%alloc(this%D, this%n_ao, this%n_ao)
!
   end subroutine initialize_density_cholesky_orbital_tool
!
!
   subroutine destruct_density_cholesky_orbital_tool(this)
!!
!!    Destruct density
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(cholesky_orbital_tool), intent(inout)  :: this
!
      call mem%dealloc(this%D, this%n_ao, this%n_ao)
!
   end subroutine destruct_density_cholesky_orbital_tool
!
!
   subroutine set_density_from_density_cholesky_orbital_tool(this, D, rank, factor)
!!
!!    Set density from density
!!    Written by Sarai D. Folkestad, 2021
!!
!!    The optional 'factor' allows for scaling. The density should be idempotent.
!!
!
      use array_initialization, only: copy_and_scale
!
      implicit none
!
      class(cholesky_orbital_tool), intent(inout)              :: this
      real(dp), dimension(this%n_ao, this%n_ao), intent(in)    :: D
!
      integer, intent(in) :: rank
!
      real(dp), intent(in), optional :: factor
!
      real(dp) :: local_factor
!
      local_factor = one
      if (present(factor)) local_factor = factor
!
      call copy_and_scale(local_factor, D, this%D, this%n_ao**2)
!
      this%initial_rank = rank
      this%rank         = rank
!
   end subroutine set_density_from_density_cholesky_orbital_tool
!
!
   subroutine set_density_from_orbitals_cholesky_orbital_tool(this, C, n)
!!
!!    Set density from orbitals
!!    Written by Sarai D. Folkestad, 2021
!!
!!    C - orbital coefficients (n_ao x n)
!!
      implicit none
!
      class(cholesky_orbital_tool), intent(inout) :: this
!
      integer, intent(in) :: n
!
      real(dp), dimension(this%n_ao, n), intent(inout) :: C
!
      call dgemm('N', 'T',    &
                  this%n_ao,  &
                  this%n_ao,  &
                  n,          &
                  one,        &
                  C,          &
                  this%n_ao,  &
                  C,          &
                  this%n_ao,  &
                  zero,       &
                  this%D,     &
                  this%n_ao)
!
      this%initial_rank = n
      this%rank         = n
!
   end subroutine set_density_from_orbitals_cholesky_orbital_tool
!
!
   subroutine restricted_decomposition_cholesky_orbital_tool(this,         &
                                                             C,            &
                                                             n_orbitals,   &
                                                             n_active_aos, &
                                                             first_active_ao)
!!
!!    Restricted decomposition
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Performs a partial and restricted Cholesky decomposition of the density.
!!    See J. Chem. Phys. 132, 204105 (2010)
!!
!!    After this routine the density becomes
!!
!!       D = D - CC^T
!!
!!    and the rank is reduced accordingly
!!
      implicit none
!
      class(cholesky_orbital_tool),                 intent(inout) :: this
      integer,                                      intent(in)    :: first_active_ao, n_active_aos
      real(dp), dimension(this%n_ao, n_active_aos), intent(out)   :: C
      integer,                                      intent(out)   :: n_orbitals
!
      integer, dimension(:), allocatable :: active_ao_list
      integer :: I
!
      call mem%alloc(active_ao_list, n_active_aos)
!
!$omp parallel do private(I)
      do I = 1, n_active_aos
!
         active_ao_list(I) = I + first_active_ao - 1
!
      enddo
!$omp end parallel do
!
      call this%cholesky_decomposition_limited_diagonal(this%D,            &
                                                         C,                &
                                                         this%n_ao,        &
                                                         n_orbitals,       &
                                                         this%threshold,   &
                                                         n_active_aos,     &
                                                         active_ao_list)
!
      call mem%dealloc(active_ao_list, n_active_aos)
!
!     Compute rank of density after decomposition
!
      this%rank = this%rank - n_orbitals
!
   end subroutine restricted_decomposition_cholesky_orbital_tool
!
!
   subroutine full_decomposition_cholesky_orbital_tool(this, C, n_orbitals)
!!
!!    Full decomposition
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Full Cholesky decomposition of the density. Terminates when all diagonals are
!!    below the threshold
!!
!!    After this routine the density should have zero rank
!
      implicit none
!
      class(cholesky_orbital_tool),                      intent(inout) :: this
      real(dp), dimension(this%n_ao, this%initial_rank), intent(out)   :: C
      integer,                                           intent(out)   :: n_orbitals
!
      integer, dimension(:), allocatable :: keep_vectors
!
      real(dp), parameter :: threshold = 1.0d-6
!
      call mem%alloc(keep_vectors, this%n_ao)
!
      call this%standard_cholesky_decomposition(this%D, C,        &
                                           this%n_ao, n_orbitals, &
                                           threshold, keep_vectors)
!
      call mem%dealloc(keep_vectors, this%n_ao)
!
!     Rank of D should now be 0
!
      this%rank = this%rank - n_orbitals
!
      if (this%rank .ne. 0) call output%error_msg('in Cholesky decomposition of the density')
!
   end subroutine full_decomposition_cholesky_orbital_tool
!
!
   subroutine cleanup_cholesky_orbital_tool(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(cholesky_orbital_tool), intent(inout) :: this
!
      call this%destruct_density()
!
      if (this%rank .lt. 0) call output%error_msg('in Cholesky decomposition of the density')
!
   end subroutine cleanup_cholesky_orbital_tool
!
!
   subroutine cholesky_decomposition_limited_diagonal(matrix, cholesky_vectors, dim_, &
                                                     n_vectors, threshold, n_included_diagonals, &
                                                     included_diagonals, n_vectors_requested)
!!
!!    Cholesky decomposition limited diagonal,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Cholesky decomposition with pivots selected from a subset of the diagonals.
!!
!!    Routine is used for decomposition of density to construct active orbitals.
!!
!!    The number of pivots may specified through the optional
!!    argument n_vectors_requested
!!
!!    On exit matrix_xy = matrix_xy - sum_J L_xJ*LyJ
!!
      implicit none
!
      integer, intent(in) :: dim_, n_included_diagonals
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in), optional :: n_vectors_requested
!
      real(dp), dimension(dim_, dim_), intent(inout) :: matrix
      real(dp), dimension(dim_, n_included_diagonals), intent(out) :: cholesky_vectors
!
      integer, dimension(n_included_diagonals), intent(in) :: included_diagonals
!
      integer, dimension(:), allocatable :: used_diag
!
      integer :: i, j, index_max, n_max_pivots
      real(dp) :: max_diagonal
!
      real(dp), dimension(:), allocatable :: diagonal
      real(dp), dimension(:), allocatable :: temp_cholesky_vector
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      n_max_pivots = n_included_diagonals
!
      if (present(n_vectors_requested)) n_max_pivots = n_vectors_requested
!
      call mem%alloc(diagonal, dim_)
!
      do i = 1, dim_
!
         diagonal(i) = matrix(i, i)
!
      enddo
!
      cholesky_vectors = zero
!
      call mem%alloc(used_diag, dim_)
      used_diag = 0
!
      do i = 1, n_max_pivots
!
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, n_included_diagonals
!
            if (abs(diagonal(included_diagonals(j))) .gt. abs(max_diagonal)) then
!
               max_diagonal = diagonal(included_diagonals(j))
               index_max    = included_diagonals(j)
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               call output%error_msg('Found negative diagonal in cholesky decomposition.')
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
!
            call mem%dealloc(used_diag, dim_)
            call mem%dealloc(diagonal, dim_)
!
!           On exit, cholesky vectors subtracted from matrix
!
            call dgemm('N', 'T',          &
                        dim_,             &
                        dim_,             &
                        n_vectors,        &
                        -one,             &
                        cholesky_vectors, &
                        dim_,             &
                        cholesky_vectors, &
                        dim_,             &
                        one,              &
                        matrix,           &
                        dim_)
!
            return
!
         else
!
            used_diag(n_vectors) = index_max
!
         endif
!
!        Cholesky vectors
!
         cholesky_vectors(:, n_vectors) = matrix(:, index_max)
!
         if (n_vectors .gt. 1) then
!
            call mem%alloc(temp_cholesky_vector, n_vectors - 1)
            temp_cholesky_vector(:) = cholesky_vectors(index_max, 1 : n_vectors - 1)
!
            call dgemm('N', 'T',                         &
                        dim_,                            &
                        1,                               &
                        n_vectors - 1,                   &
                        -one,                            &
                        cholesky_vectors,                &
                        dim_,                            &
                        temp_cholesky_vector,            &
                        1,                               &
                        one,                             &
                        cholesky_vectors(1, n_vectors),  &
                        dim_)
!
            call mem%dealloc(temp_cholesky_vector, n_vectors - 1)
!
         endif
!
         do j = 1, n_vectors - 1
!
            cholesky_vectors(used_diag(j), n_vectors) = zero
!
         enddo
!
         call dscal(dim_, one/sqrt(max_diagonal), cholesky_vectors(1, n_vectors), 1)
!
         do j = 1, dim_
!
            diagonal(j) = diagonal(j) - cholesky_vectors(j, n_vectors)**2
!
         enddo
!
         diagonal(index_max) = zero
!
      enddo
!
      call mem%dealloc(used_diag, dim_)
      call mem%dealloc(diagonal, dim_)
!
!     On exit, cholesky vectors subtracted from matrix
!
      call dgemm('N', 'T',          &
                  dim_,             &
                  dim_,             &
                  n_vectors,        &
                  -one,             &
                  cholesky_vectors, &
                  dim_,             &
                  cholesky_vectors, &
                  dim_,             &
                  one,              &
                  matrix,           &
                  dim_)
!
   end subroutine cholesky_decomposition_limited_diagonal
!
!
   subroutine standard_cholesky_decomposition(matrix, cholesky_vectors, dim_, &
                                          n_vectors, threshold, used_diag)
!!
!!    Standard Cholesky decomposition,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Partial pivoting (column) Cholesky decomposition algorithm.
!!    Selects largest diagonal element as pivot.
!!
!!    Completes when all diagonals are below threshold.
!!
      implicit none
!
      integer, intent(in) :: dim_
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(dim_, dim_), intent(inout) :: matrix
      real(dp), dimension(dim_, dim_), intent(out) :: cholesky_vectors
!
      integer, dimension(dim_), intent(out) :: used_diag
!
      integer :: i, j, index_max
      real(dp) :: max_diagonal, min_diagonal
!
      real(dp), dimension(:), allocatable :: diagonal
      real(dp), dimension(:), allocatable :: temp_cholesky_vector
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      used_diag = 0
!
!     Looping over the number of cholesky vectors
!
      call mem%alloc(diagonal, dim_)
!
      do i = 1, dim_
!
         diagonal(i) = matrix(i, i)
!
      enddo
!
      do i = 1, dim_
!
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, dim_
!
            if (abs(diagonal(j)) .gt. abs(max_diagonal)) then
!
               max_diagonal = diagonal(j)
               index_max    = j
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               call output%error_msg('Found negative diagonal in cholesky decomposition.')
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
!
            min_diagonal = 1.0D10
!
            do j = 1, dim_
!
            call dgemm('N', 'T',          &
                        dim_,             &
                        dim_,             &
                        n_vectors,        &
                        -one,             &
                        cholesky_vectors, &
                        dim_,             &
                        cholesky_vectors, &
                        dim_,             &
                        one,              &
                        matrix,           &
                        dim_)
               if (diagonal(j) .lt. min_diagonal) min_diagonal = diagonal(j)

!
            enddo
!
            call output%printf('n', 'The smallest diagonal after decomposition &
                               &is: (e12.4)', reals=[min_diagonal], fs='(/t6,a)')
            call mem%dealloc(diagonal, dim_)
!
            return
!
         else
!
            used_diag(n_vectors) = index_max
!
         endif
!
!        Cholesky vectors
!
         cholesky_vectors(:, n_vectors) = matrix(:, index_max)
!
         if (n_vectors .gt. 1) then
!
            call mem%alloc(temp_cholesky_vector, n_vectors - 1)
            temp_cholesky_vector(:) = cholesky_vectors(index_max, 1 : n_vectors - 1)
!
            call dgemm('N', 'T',                         &
                        dim_,                            &
                        1,                               &
                        n_vectors - 1,                   &
                        -one,                            &
                        cholesky_vectors,                &
                        dim_,                            &
                        temp_cholesky_vector,            &
                        1,                               &
                        one,                             &
                        cholesky_vectors(1, n_vectors),  &
                        dim_)
!
            call mem%dealloc(temp_cholesky_vector, n_vectors - 1)
!
         endif
!
         do j = 1, n_vectors - 1
!
            cholesky_vectors(used_diag(j), n_vectors) = zero
!
         enddo
!
         call dscal(dim_, one/sqrt(max_diagonal), cholesky_vectors(1, n_vectors), 1)
!
         do j = 1, dim_
!
            diagonal(j) = diagonal(j) - cholesky_vectors(j, n_vectors)**2
!
         enddo
!
         diagonal(index_max) = zero
!
         do j = 1, dim_
!
            matrix(j,index_max) = 0.0D0
            matrix(index_max,j) = 0.0D0
!
         enddo
!
      enddo
!
      min_diagonal = 1.0D10
!
      do j = 1, dim_
!
            if (diagonal(j) .lt. min_diagonal) min_diagonal = diagonal(j)

!
      enddo
!
      call output%printf('n', 'The smallest diagonal after decomposition is: (e12.4)', &
                         reals=[min_diagonal], fs='(/t6,a)')
!
      call mem%dealloc(diagonal, dim_)
!
   end subroutine standard_cholesky_decomposition
!
!
end module cholesky_orbital_tool_class
