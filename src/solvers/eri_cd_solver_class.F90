module eri_cd_solver_class
!
!!
!!    Electronic repulsion ao_integrals solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
!
   use index
   use reordering
   use array_utilities
   use array_analysis
!
   use file_class
   use disk_manager_class
!
   use io_utilities
!
   use cholesky_array_list_class
!
   use wavefunction_class
!
   implicit none
!
   type :: eri_cd_solver
!
      real(dp) :: threshold   = 1.0D-8
      real(dp) :: span        = 1.0D-2
!
      integer(i15) :: max_qual = 1000
!
      integer(i15) :: iteration = 0
!
      logical :: one_center         = .false.
      logical :: construct_vectors  = .true.
!
      type(file) :: diagonal_info_target, diagonal_info_one_center
      type(file) :: cholesky_aux, cholesky_aux_inverse, cholesky_ao_vectors, cholesky_ao_vectors_info
      type(file) :: basis_shell_data
!
      integer(i15) :: n_cholesky
      integer(i15) :: n_s, n_sp, n_ao, n_aop
!
   contains
!
      procedure :: initialize => initialize_eri_cd_solver
      procedure :: solve      => solve_eri_cd_solver
      procedure :: finalize   => finalize_eri_cd_solver
!
      procedure :: invert_overlap_cholesky_vecs                => invert_overlap_cholesky_vecs_eri_cd_solver
      procedure :: cholesky_vecs_diagonal_test                 => cholesky_vecs_diagonal_test_eri_cd_solver
      procedure :: construct_significant_diagonal              => construct_significant_diagonal_eri_cd_solver
      procedure :: construct_significant_diagonal_atomic       => construct_significant_diagonal_atomic_eri_cd_solver
      procedure :: determine_auxilliary_cholesky_basis         => determine_auxilliary_cholesky_basis_eri_cd_solver
      procedure :: construct_overlap_cholesky_vecs             => construct_overlap_cholesky_vecs_eri_cd_solver
      procedure :: construct_cholesky_vectors                  => construct_cholesky_vectors_eri_cd_solver
!
      procedure :: read_info  => read_info_eri_cd_solver
!
   end type eri_cd_solver
!
!
contains
!
!
   subroutine initialize_eri_cd_solver(solver, system)
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      type(molecular_system) :: system
!
      if (requested_section('cholesky')) then
         call solver%read_info()
      endif
!
      solver%n_aop   = system%get_n_aos()*(system%get_n_aos()+1)/2 ! Number of ao pairs packed
      solver%n_ao    = system%get_n_aos()
      solver%n_s     = system%get_n_shells()
      solver%n_sp    = solver%n_s*(solver%n_s + 1)/2         ! Number of shell pairs packed
!
!
      call solver%diagonal_info_target%init('target_diagonal', 'sequential', 'unformatted')
      call solver%cholesky_aux%init('cholesky_aux', 'sequential', 'unformatted')
      call solver%cholesky_aux_inverse%init('cholesky_aux_inverse', 'sequential', 'unformatted')
      call solver%basis_shell_data%init('basis_shell_info', 'sequential', 'unformatted')
!
      if (solver%one_center) then
!
         call solver%diagonal_info_one_center%init('one_center_diagonal', 'sequential', 'unformatted')
!
      endif
!
        call initialize_coulomb()
        call initialize_kinetic()
        call initialize_nuclear()
        call initialize_overlap()
!
   end subroutine initialize_eri_cd_solver
!
!
   subroutine solve_eri_cd_solver(solver, system, screening_vector)
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system) :: system
!
      real(dp), dimension(solver%n_ao,1), optional :: screening_vector
!
      real(dp):: s_determine_basis, e_determine_basis, s_build_vectors, e_build_vectors, omp_get_wtime
!
      write(output%unit, '(/a/)') ':: Cholesky decomposition of two-electron ao integrals'
      flush(output%unit)
!
      write(output%unit, '(a20, i10)')'Number of aos:      ', solver%n_ao
      write(output%unit, '(a20, i10)')'Number of ao pairs: ', solver%n_aop
      write(output%unit, '(a20, i10)')'Number of shells:   ', solver%n_s
!
      write(output%unit, '(/a21, e12.4)') 'Target threshold is: ', solver%threshold
      write(output%unit, '(a21, e12.4/)') 'Span factor:         ', solver%span
      flush(output%unit)
!
      s_determine_basis = omp_get_wtime()
!
      if (present(screening_vector)) then
!
         call solver%construct_significant_diagonal(system, screening_vector)
!
      else 
!
         call solver%construct_significant_diagonal(system)
!
      endif
!
      if (solver%one_center) then
!
         if (present(screening_vector)) then
!
            call solver%construct_significant_diagonal_atomic(system, screening_vector)
!
         else 
!
            call solver%construct_significant_diagonal_atomic(system)
!
         endif
!
         call solver%determine_auxilliary_cholesky_basis(system, solver%diagonal_info_one_center)
!
      else
!
         call solver%determine_auxilliary_cholesky_basis(system, solver%diagonal_info_target)
!
      endif
!
      call solver%construct_overlap_cholesky_vecs(system)
      call solver%invert_overlap_cholesky_vecs()
!
      e_determine_basis = omp_get_wtime()
      write(output%unit, '(/a30, f11.2)')'Wall time to determine basis: ', &
                                  e_determine_basis - s_determine_basis
!
      if (solver%construct_vectors) then
!
         s_build_vectors = omp_get_wtime()
!
         call solver%construct_cholesky_vectors(system)
!
         e_build_vectors = omp_get_wtime()
         write(output%unit, '(/a37, f11.2)')'Wall time to build vectors and test: ', &
                                  e_build_vectors - s_build_vectors
!
      endif
!
   end subroutine solve_eri_cd_solver
!
!
   subroutine finalize_eri_cd_solver(solver)
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
   end subroutine finalize_eri_cd_solver
!
!
   subroutine construct_significant_diagonal_eri_cd_solver(solver, system, screening_vector)
!!
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      class(molecular_system) :: system
!
      real(dp), dimension(solver%n_ao,1), optional :: screening_vector
!
      integer(i15) ::sp, n_sig_aop, n_sig_sp, current_sig_sp
      integer(i15), dimension(:,:), allocatable :: sp_index, sig_sp_index, ao_offsets
!
      real(dp), dimension(:,:), allocatable :: g_AB_AB, D_AB, D_AB_screen, D_xy, screening_vector_local, screening_vector_reduced
!
      integer(i15) :: x, y, xy, xy_packed, first_sig_sp, first_sig_aop, A, B, I
!
      type(interval) :: A_interval, B_interval
!
      logical, dimension(:), allocatable :: sig_sp
!
      call mem%alloc(screening_vector_local, solver%n_aop, 1)
!
      if (present(screening_vector)) then
!
         screening_vector_local = screening_vector
!
      else
!
         screening_vector_local = one
!
      endif
!
!     Prepare for pre-screening
!
      call mem%alloc_int(sp_index, solver%n_sp, 2)
!
      sp = 0        ! Shell pair number
!
      do B = 1, solver%n_s
         do A = B, solver%n_s
!
            sp = sp + 1
!
            sp_index(sp, 1) = A
            sp_index(sp, 2) = B
!
         enddo
      enddo
!
!     Pre-screening of full diagonal
!
      allocate(sig_sp(solver%n_sp))
      sig_sp = .false.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, g_AB_AB, D_AB, D_AB_screen) &
!$omp shared(sig_sp) &
!$omp schedule(guided)
      do I = 1, solver%n_sp
!
         A = sp_index(I, 1)
         B = sp_index(I, 2)
!
         A_interval = system%shell_limits(A)
         B_interval = system%shell_limits(B)
!
!        Construct diagonal D_AB for the given shell pair
!
         call mem%alloc(g_AB_AB, &
                  (A_interval%size)*(B_interval%size), &
                  (A_interval%size)*(B_interval%size))
!
         g_AB_AB = zero
         call system%ao_integrals%get_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
         call mem%alloc(D_AB, (A_interval%size)*(B_interval%size), 1)
         call mem%alloc(D_AB_screen, (A_interval%size)*(B_interval%size), 1)
!
         do x = 1, (A_interval%size)
            do y = 1, (B_interval%size)
!
               xy = (A_interval%size)*(y-1)+x
!
               D_AB_screen(xy, 1) = g_AB_AB(xy, xy)&
                           *screening_vector_local(x + A_interval%first - 1,1)&
                           *screening_vector_local(y + B_interval%first - 1,1)
               D_AB(xy, 1) = g_AB_AB(xy, xy)
!
            enddo
         enddo
!
         call mem%dealloc(g_AB_AB, &
                  (A_interval%size)*(B_interval%size), &
                  (A_interval%size)*(B_interval%size))
!
!        Determine whether shell pair is significant
!
         sig_sp(I) = (is_significant(D_AB, (A_interval%size)*(B_interval%size), solver%threshold) .and. &
                      is_significant(D_AB_screen, (A_interval%size)*(B_interval%size), solver%threshold))
!
         call mem%dealloc(D_AB, (A_interval%size)*(B_interval%size), 1)
         call mem%dealloc(D_AB_screen, (A_interval%size)*(B_interval%size), 1)
!
      enddo
!$omp end parallel do
!
      n_sig_aop = 0 ! Number of significant AO pairs
      n_sig_sp  = 0 ! Number of significant shell pairs
!
      do I = 1, solver%n_sp
!
         if (sig_sp(I)) then
!
            A = sp_index(I, 1)
            B = sp_index(I, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            n_sig_aop = n_sig_aop + &
                           get_size_sp(A_interval, B_interval)
!
            n_sig_sp = n_sig_sp + 1
!
         endif
!
      enddo
!
      write(output%unit, '(/a)')'Reduction of shell pairs:'
      write(output%unit, '(a33, 2x, i9)')'Total number of shell pairs:     ', solver%n_sp
      write(output%unit, '(a33, 2x, i9)')'Significant shell pairs:         ', n_sig_sp
      write(output%unit, '(a33, 2x, i9)')'Significant ao pairs:            ', n_sig_aop
      flush(output%unit)
!
!     Prepare for construction of diagonal and screening vector
!     
      call mem%alloc_int(ao_offsets, n_sig_sp, 1)
      ao_offsets = 0
!
      current_sig_sp = 0
!
      call mem%alloc_int(sig_sp_index, n_sig_sp, 2)
!
      do I = 1, solver%n_sp
!
         if (sig_sp(I)) then
!
            current_sig_sp = current_sig_sp + 1
!
            A = sp_index(I, 1)
            B = sp_index(I, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            sig_sp_index(current_sig_sp, 1) = A
            sig_sp_index(current_sig_sp, 2) = B
!
            if (current_sig_sp .lt. n_sig_sp) then
!
               ao_offsets(current_sig_sp + 1, 1) = ao_offsets(current_sig_sp, 1) + &
                           get_size_sp(A_interval, B_interval)
!
            endif
!
         endif
!
      enddo
!
      call mem%dealloc_int(sp_index, solver%n_sp, 2)
!
!     Construct significant diagonal and screening vector
!
      call mem%alloc(D_xy, n_sig_aop, 1)
      D_xy = zero
!
      call mem%alloc(screening_vector_reduced, n_sig_aop, 1)
      screening_vector_reduced = zero
!
!     Note: allocated with length n_significant_sp + 1, last element is used for n_significant_aop
!     This is convenient because significant_sp_to_first_significant_aop will be used to calculate lengths.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, xy_packed, g_AB_AB) &
!$omp shared(D_xy, screening_vector_reduced, ao_offsets) &
!$omp schedule(guided)
      do I = 1, n_sig_sp
!
         A = sig_sp_index(I, 1)
         B = sig_sp_index(I, 2)
!
         A_interval = system%shell_limits(A)
         B_interval = system%shell_limits(B)
!
         call mem%alloc(g_AB_AB, &
               (A_interval%size)*(B_interval%size), &
               (A_interval%size)*(B_interval%size))
!
         g_AB_AB = zero
         call system%ao_integrals%get_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
         if (A .eq. B) then
!
            do x = 1, A_interval%size
               do y = x, B_interval%size
!
                  xy = A_interval%size*(y - 1) + x
                  xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                  D_xy(xy_packed + ao_offsets(I, 1), 1) = g_AB_AB(xy, xy)
                  screening_vector_reduced(xy_packed + ao_offsets(I, 1), 1) = &
                                          screening_vector_local(x + A_interval%first - 1, 1)*&
                                          screening_vector_local(y + B_interval%first - 1, 1)
!
               enddo
            enddo
!
         else ! A ≠ B
!
            do x = 1, (A_interval%size)
               do y = 1, (B_interval%size)
!
                  xy = A_interval%size*(y - 1) + x
                  D_xy(xy + ao_offsets(I, 1), 1) = g_AB_AB(xy,xy)
                  screening_vector_reduced(xy + ao_offsets(I, 1), 1) = &
                                      screening_vector_local(x + A_interval%first - 1, 1)*&
                                      screening_vector_local(y + B_interval%first - 1, 1)
!
               enddo
            enddo
!
         endif
!
         call mem%dealloc(g_AB_AB, &
               (A_interval%size)*(B_interval%size), &
               (A_interval%size)*(B_interval%size))
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(screening_vector_local, solver%n_aop, 1)
      call mem%dealloc_int(sig_sp_index, n_sig_sp, 2)
!
!
!     Write screening_info_file containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call disk%open_file(solver%diagonal_info_target, 'write')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
      write(solver%diagonal_info_target%unit) sig_sp
      write(solver%diagonal_info_target%unit) D_xy
      write(solver%diagonal_info_target%unit) screening_vector_reduced
!
      call disk%close_file(solver%diagonal_info_target)
!
      deallocate(sig_sp)
      call mem%dealloc(D_xy, n_sig_aop, 1)
      call mem%dealloc(screening_vector_reduced, n_sig_aop, 1)
!
   end subroutine construct_significant_diagonal_eri_cd_solver
!
!
   subroutine construct_significant_diagonal_atomic_eri_cd_solver(solver, system, screening_vector)
!!
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      class(molecular_system) :: system
!
      real(dp), dimension(solver%n_ao,1), optional :: screening_vector
!
      integer(i15) ::sp, n_sig_aop, n_sig_sp, current_sig_sp
      integer(i15), dimension(:,:), allocatable :: sp_index, sig_sp_index, ao_offsets
!
      real(dp), dimension(:,:), allocatable :: g_AB_AB, D_AB, D_AB_screen, D_xy, screening_vector_local, screening_vector_reduced
!
      integer(i15) :: x, y, xy, xy_packed, first_sig_sp, first_sig_aop, A, B, I
!
      type(interval) :: A_interval, B_interval
!
      logical, dimension(:), allocatable :: sig_sp
!
      call mem%alloc(screening_vector_local, solver%n_aop, 1)
!
      if (present(screening_vector)) then
!
         screening_vector_local = screening_vector
!
      else
!
         screening_vector_local = one
!
      endif
!
!     Prepare for pre-screening
!
      call mem%alloc_int(sp_index, solver%n_sp, 2)
!
      sp = 0        ! Shell pair number
!
      do B = 1, solver%n_s
         do A = B, solver%n_s
!
            sp = sp + 1
!
            sp_index(sp, 1) = A
            sp_index(sp, 2) = B
!
         enddo
      enddo
!
!     Pre-screening of full diagonal
!
      allocate(sig_sp(solver%n_sp))
      sig_sp = .false.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, g_AB_AB, D_AB, D_AB_screen) &
!$omp shared(sig_sp) &
!$omp schedule(guided)
      do I = 1, solver%n_sp
!
         A = sp_index(I, 1)
         B = sp_index(I, 2)
!
         if (system%shell_to_atom(B) == system%shell_to_atom(A)) then
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
!           Construct diagonal D_AB for the given shell pair
!
            call mem%alloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
            g_AB_AB = zero
            call system%ao_integrals%get_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
            call mem%alloc(D_AB, (A_interval%size)*(B_interval%size), 1)
            call mem%alloc(D_AB_screen, (A_interval%size)*(B_interval%size), 1)
!
            do x = 1, (A_interval%size)
               do y = 1, (B_interval%size)
!
                  xy = (A_interval%size)*(y-1)+x
!
                  D_AB_screen(xy, 1) = g_AB_AB(xy, xy)&
                              *screening_vector_local(x + A_interval%first - 1,1)&
                              *screening_vector_local(y + B_interval%first - 1,1)
!  
                  D_AB(xy, 1) = g_AB_AB(xy, xy)
!
               enddo
            enddo
!
            call mem%dealloc(g_AB_AB, &
                  (A_interval%size)*(B_interval%size), &
                  (A_interval%size)*(B_interval%size))
!
!           Determine whether shell pair is significant
!
            sig_sp(I) = (is_significant(D_AB, (A_interval%size)*(B_interval%size), solver%threshold) .and. &
                         is_significant(D_AB_screen, (A_interval%size)*(B_interval%size), solver%threshold))
!
            call mem%dealloc(D_AB, (A_interval%size)*(B_interval%size), 1)
            call mem%dealloc(D_AB_screen, (A_interval%size)*(B_interval%size), 1)
!
         endif
!
      enddo
!$omp end parallel do
!
      n_sig_aop = 0 ! Number of significant AO pairs
      n_sig_sp  = 0 ! Number of significant shell pairs
!
      do I = 1, solver%n_sp
!
         if (sig_sp(I)) then
!
            A = sp_index(I, 1)
            B = sp_index(I, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            n_sig_aop = n_sig_aop + &
                           get_size_sp(A_interval, B_interval)
!
            n_sig_sp = n_sig_sp + 1
!
         endif
!
      enddo
!
      write(output%unit, '(/a)')'Reduction of shell pairs:'
      write(output%unit, '(a33, 2x, i9)')'Total number of shell pairs:     ', solver%n_sp
      write(output%unit, '(a33, 2x, i9)')'Significant shell pairs:         ', n_sig_sp
      write(output%unit, '(a33, 2x, i9)')'Significant ao pairs:            ', n_sig_aop
      flush(output%unit)
!
!     Prepare for construction of diagonal and screening vector
!     
      call mem%alloc_int(ao_offsets, n_sig_sp, 1)
      ao_offsets = 0
!
      current_sig_sp = 0
!
      call mem%alloc_int(sig_sp_index, n_sig_sp, 2)
!
      do I = 1, solver%n_sp
!
         if (sig_sp(I)) then
!
            current_sig_sp = current_sig_sp + 1
!
            A = sp_index(I, 1)
            B = sp_index(I, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            sig_sp_index(current_sig_sp, 1) = A
            sig_sp_index(current_sig_sp, 2) = B
!
            if (current_sig_sp .lt. n_sig_sp) then
!
               ao_offsets(current_sig_sp + 1, 1) = ao_offsets(current_sig_sp, 1) + &
                           get_size_sp(A_interval, B_interval)
!
            endif
!
         endif
!
      enddo
!
      call mem%dealloc_int(sp_index, solver%n_sp, 2)
!
!     Construct significant diagonal and screening vector
!
      call mem%alloc(D_xy, n_sig_aop, 1)
      D_xy = zero
!
      call mem%alloc(screening_vector_reduced, n_sig_aop, 1)
      screening_vector_reduced = zero
!
!     Note: allocated with length n_significant_sp + 1, last element is used for n_significant_aop
!     This is convenient because significant_sp_to_first_significant_aop will be used to calculate lengths.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, xy_packed, g_AB_AB) &
!$omp shared(D_xy, screening_vector_reduced, ao_offsets) &
!$omp schedule(guided)
      do I = 1, n_sig_sp
!
         A = sig_sp_index(I, 1)
         B = sig_sp_index(I, 2)
!
         A_interval = system%shell_limits(A)
         B_interval = system%shell_limits(B)
!
         call mem%alloc(g_AB_AB, &
               (A_interval%size)*(B_interval%size), &
               (A_interval%size)*(B_interval%size))
!
         g_AB_AB = zero
         call system%ao_integrals%get_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
         if (A .eq. B) then
!
            do x = 1, A_interval%size
               do y = x, B_interval%size
!
                  xy = A_interval%size*(y - 1) + x
                  xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                  D_xy(xy_packed + ao_offsets(I, 1), 1) = g_AB_AB(xy, xy)
                  screening_vector_reduced(xy_packed + ao_offsets(I, 1), 1) = &
                                          screening_vector_local(x + A_interval%first - 1, 1)*&
                                          screening_vector_local(y + B_interval%first - 1, 1)
!
               enddo
            enddo
!
         else ! A ≠ B
!
            do x = 1, (A_interval%size)
               do y = 1, (B_interval%size)
!
                  xy = A_interval%size*(y - 1) + x
                  D_xy(xy + ao_offsets(I, 1), 1) = g_AB_AB(xy,xy)
                  screening_vector_reduced(xy + ao_offsets(I, 1), 1) = &
                                      screening_vector_local(x + A_interval%first - 1, 1)*&
                                      screening_vector_local(y + B_interval%first - 1, 1)
!
               enddo
            enddo
!
         endif
!
         call mem%dealloc(g_AB_AB, &
               (A_interval%size)*(B_interval%size), &
               (A_interval%size)*(B_interval%size))
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(screening_vector_local, solver%n_aop, 1)
      call mem%dealloc_int(sig_sp_index, n_sig_sp, 2)
!
!
!     Write screening_info_file containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call disk%open_file(solver%diagonal_info_target, 'write')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
      write(solver%diagonal_info_target%unit) sig_sp
      write(solver%diagonal_info_target%unit) D_xy
      write(solver%diagonal_info_target%unit) screening_vector_reduced
!
      call disk%close_file(solver%diagonal_info_target)
!
      deallocate(sig_sp)
      call mem%dealloc(D_xy, n_sig_aop, 1)
      call mem%dealloc(screening_vector_reduced, n_sig_aop, 1)
!
   end subroutine construct_significant_diagonal_atomic_eri_cd_solver
!
!
   subroutine determine_auxilliary_cholesky_basis_eri_cd_solver(solver, system, diagonal_info)
!!
!!    ....
!!
!!
      implicit none
!
      class(eri_cd_solver), intent(inout) :: solver
!
      class(molecular_system), intent(in) :: system
!
      type(file), intent(in) :: diagonal_info
!
!     Local variables
!
!     Integers
      integer(i15) :: n_sig_sp, n_sig_aop
      integer(i15) :: n_new_sig_sp, n_new_sig_aop, current_new_sig_sp
      integer(i15) :: n_qual_sp, n_qual_aop, n_previous_qual_aop, n_qual_aop_in_sp
      integer(i15) :: sp, current_sig_sp
      integer(i15) :: first_sig_aop, last_sig_aop, aop
      integer(i15) :: A, B, AB, AB_sp
      integer(i15) :: C, D, CD_sp
      integer(i15) :: I, J
      integer(i15) :: w, x, y, z
      integer(i15) :: xy, xy_packed, xy_max, wx, wx_packed, yz
      integer(i15) :: sig_neg
      integer(i15) :: first, last
      integer(i15) :: first_x, first_y
      integer(i15) :: current_qual, qual
      integer(i15) :: n_cholesky_in_node, n_new_cholesky
      integer(i15) :: n_sp_in_basis
      integer(i15) :: sig_sp_counter
      integer(i15) :: sp_in_basis
      integer(kind=4) :: omp_get_thread_num
!
!     Integer allocatable arrays
!
      integer(i15), dimension(:,:), allocatable :: sig_sp_to_first_sig_aop       ! Maps significant shell pair to first ao pair
      integer(i15), dimension(:,:), allocatable :: new_sig_sp_to_first_sig_aop   ! Maps significant shell pair to first ao pair
      integer(i15), dimension(:,:), allocatable :: sig_sp_to_shells              ! Maps significant shell pair to shells
      integer(i15), dimension(:,:), allocatable :: new_sig_sp_to_shells          ! Maps significant shell pair to shells
      integer(i15), dimension(:,:), allocatable :: sig_aop_to_aos                ! Maps significant ao pair to aos
      integer(i15), dimension(:,:), allocatable :: new_sig_aop_to_aos            ! Maps significant ao pair to aos
      integer(i15), dimension(:,:), allocatable :: sorted_max_sig_sp             ! Index array for sorting shell pairs according to their maximum values
      integer(i15), dimension(:,:), allocatable :: qual_sp                       ! List of qualified shell pairs
      integer(i15), dimension(:,:), allocatable :: qual_sp_copy                  ! List of qualified shell pairs, copy used to reduce size
      integer(i15), dimension(:,:), allocatable :: qual_aop                      ! List of qualified ao pairs
      integer(i15), dimension(:,:), allocatable :: qual_aop_copy                 ! List of qualified ao pairs, copy used to reduce size
      integer(i15), dimension(:,:), allocatable :: sorted_qual_aop_in_sp_indices ! Index array for sorting the qualified ao pairs in shell pair
      integer(i15), dimension(:,:), allocatable :: n_qual_aop_in_prev_sps        ! Offsets for omp-loop, number of qualified ao pairs in preceding shell pair
      integer(i15), dimension(:,:), allocatable :: cholesky_basis                ! ao and ao pair indices of the elements of the cholesky basis
      integer(i15), dimension(:,:), allocatable :: cholesky_basis_new            ! ao and ao pair indices of the elements of the cholesky basis, written to file at end of routine
      integer(i15), dimension(:,:), allocatable :: qual_max                      ! Index list containing order in which qualified diagonals are selected in decomposition
      integer(i15), dimension(:,:), allocatable :: sig_sp_to_previous_sig_sp     ! Maps significant shell pair indices to significant shell pair indices of last iteration, used for reduction
      integer(i15), dimension(:,:), allocatable :: basis_shell_info_full         ! Info on shells containing elements of the basis
      integer(i15), dimension(:,:), allocatable :: basis_shell_info              ! Info on shells containing elements of the basis, written to file at end of routine
!
!     Logicals
!
      logical :: done, write_warning, construct_more_choleskys, found
!
!     Logical allocatable arrays
!
      logical, dimension(:), allocatable :: sig_sp
      logical, dimension(:), allocatable :: new_sig_sp
!
!     Reals
!
      real(dp) :: D_max_full, D_max
!
!     Reals for timings
!
      real(dp) :: s_select_basis_time, e_select_basis_time
      real(dp) :: full_reduce_time, s_reduce_time, e_reduce_time
      real(dp) :: full_construct_time, s_construct_time, e_construct_time
!
!     Real allocatable arrays
!
      real(dp), dimension(:,:), allocatable :: D_xy                           ! Array for eri diagonal elements
      real(dp), dimension(:,:), allocatable :: D_xy_new                       ! Array for eri diagonal elements, used for reduction
      real(dp), dimension(:,:), allocatable :: approx_diagonal_accumulative   ! Array for accumulating approximate diagonal
      real(dp), dimension(:,:), allocatable :: g_wxyz                         ! Array for eri
      real(dp), dimension(:,:), allocatable :: g_AB_CD                        ! Array for eri for shell pairs AB and CD
      real(dp), dimension(:,:), allocatable :: cholesky_tmp                   ! Array used for dgemm, reordered copy of cholesky vectors of current batch of qualified
      real(dp), dimension(:,:), allocatable :: max_in_sig_sp                  ! Maximum in each significant shell pair
      real(dp), dimension(:,:), allocatable :: sorted_qual_aop_in_sp          ! Sorted qualified ao pair in shell pair
      real(dp), dimension(:,:), allocatable :: screening_vector               ! Screening vector for diagonal
      real(dp), dimension(:,:), allocatable :: screening_vector_new           ! Screening vector for diagonal, used for reduction
!
!     Real pointers
!
      real(dp), dimension(:,:), pointer     :: cholesky
      real(dp), dimension(:,:), pointer     :: cholesky_new
!
!     Intervals
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
!     Cholesky linked list
!
      type(cholesky_array_list) :: cholesky_array
!
      call cpu_time(s_select_basis_time)

!     Read diagonal info file containing (name given as argument)
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call disk%open_file(diagonal_info, 'read')
!
      rewind(diagonal_info%unit)
!
      read(diagonal_info%unit) n_sig_sp, n_sig_aop
!
      call mem%alloc(D_xy, n_sig_aop, 1)
      call mem%alloc(screening_vector, n_sig_aop, 1)
      allocate(sig_sp(solver%n_sp))
!
      read(diagonal_info%unit) sig_sp
      read(diagonal_info%unit) D_xy
      read(diagonal_info%unit) screening_vector
!
      call disk%close_file(diagonal_info)
!
!     Construct info arrays
!
      call mem%alloc_int(sig_sp_to_first_sig_aop, n_sig_sp + 1, 1)
      sig_sp_to_first_sig_aop = 0
!
      sig_sp_to_first_sig_aop(n_sig_sp + 1, 1) = n_sig_aop + 1
!
      call mem%alloc_int(sig_sp_to_shells, n_sig_sp, 2) ! [A, B]
      sig_sp_to_shells = 0
!
      call mem%alloc_int(sig_aop_to_aos, n_sig_aop, 2) ! [alpha, beta]
      sig_aop_to_aos = 0
!
!     Note: allocated with length n_significant_sp + 1, last element is used for n_significant_aop
!     This is convenient because significant_sp_to_first_significant_aop will be used to calculate lengths.
!
      sp              = 1
      current_sig_sp  = 1
      first_sig_aop   = 1
!
      do B = 1, solver%n_s
         do A = B, solver%n_s
!
            if (sig_sp(sp)) then
!
               sig_sp_to_first_sig_aop(current_sig_sp, 1) = first_sig_aop
!
               A_interval = system%shell_limits(A)
               B_interval = system%shell_limits(B)
!
               sig_sp_to_shells(current_sig_sp, 1) = A
               sig_sp_to_shells(current_sig_sp, 2) = B
!
               if (A .eq. B) then
!
                  do x = 1, A_interval%size
                     do y = 1, B_interval%size
!
                        xy = A_interval%size*(y - 1) + x
                        xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                        sig_aop_to_aos(xy_packed + first_sig_aop - 1, 1) = &
                                                      A_interval%first + x - 1
!
                        sig_aop_to_aos(xy_packed + first_sig_aop - 1, 2) = &
                                                      B_interval%first + y - 1
!
                     enddo
                  enddo
!
               else ! A ≠ B
!
                  do x = 1, (A_interval%size)
                     do y = 1, (B_interval%size)
!
                        xy = A_interval%size*(y - 1) + x
!
                        sig_aop_to_aos(xy + first_sig_aop - 1, 1) = A_interval%first + x - 1
                        sig_aop_to_aos(xy + first_sig_aop - 1, 2) = B_interval%first + y - 1
!
                     enddo
                  enddo
!
               endif
!
               first_sig_aop = first_sig_aop + get_size_sp(A_interval, B_interval)
!
               current_sig_sp = current_sig_sp + 1
!
            endif ! End of if (significant)
!
            sp = sp + 1
!
         enddo
      enddo
!
!     Determining the basis
!
      write(output%unit, '(/a)') ' - Determinig the elements of the basis'
      flush(output%unit)
!
      write(output%unit, '(/a)')&
      'Iter.  #Sign. ao pairs / shell pairs   Max diagonal    #Qualified    #Cholesky    Cholesky array size'
      write(output%unit, '(a)') &
      '-------------------------------------------------------------------------------------------------------'
      flush(output%unit)
!
      solver%iteration   = 0
      solver%n_cholesky  = 0
!
      sig_neg = 0 ! Counts number of significant neative diagonals (absolute value > 10^-10)
!
      full_reduce_time     = zero
      full_construct_time  = zero
!
      done = .false.
!
      call cholesky_array%initialize()
!
      do while (.not. done)
!
         write_warning = .true. ! Logical used to ensure warning of significant negative diagonal only appears once per iteration
!
         solver%iteration = solver%iteration + 1
!
!        Shell maximums and shell maximums indices vectors
!
         call mem%alloc(max_in_sig_sp, n_sig_sp, 1)
!
         max_in_sig_sp = zero
!
         do sp = 1, n_sig_sp
!
!           Get first and last indices of shell pair
!
            first = sig_sp_to_first_sig_aop(sp, 1)
            last  = sig_sp_to_first_sig_aop(sp + 1, 1) - 1
!
!           Determine the largest elements
!
            do I = first, last
!
               if (D_xy(I, 1) .gt. max_in_sig_sp(sp, 1)) then
!
                  max_in_sig_sp(sp, 1) = D_xy(I, 1)
!
               endif
!
            enddo
!
         enddo
!
!        Sort from largest to smallest and determine an index array of sorting
!
         call mem%alloc_int(sorted_max_sig_sp, n_sig_sp,1)
         sorted_max_sig_sp = 0
!
         call quicksort_with_index(max_in_sig_sp, sorted_max_sig_sp, n_sig_sp)
!
         D_max_full  = max_in_sig_sp(1, 1)
         n_qual_aop  = 0
         n_qual_sp   = 0
!
         call mem%dealloc(max_in_sig_sp, n_sig_sp, 1)
!
         call mem%alloc_int(qual_aop, solver%max_qual, 3)
         call mem%alloc_int(qual_sp, solver%n_sp, 3)
         qual_sp = 0
         qual_aop = 0
!
         do sp = 1, n_sig_sp
!
            current_sig_sp = sorted_max_sig_sp(sp, 1)
!
            first_sig_aop = sig_sp_to_first_sig_aop(current_sig_sp, 1)
            last_sig_aop  = sig_sp_to_first_sig_aop(current_sig_sp + 1, 1) - 1
!
            n_qual_aop_in_sp = 0
!
            do aop = first_sig_aop, last_sig_aop
!
               if ((D_xy(aop, 1) .ge. solver%span*D_max_full) .and. (n_qual_aop .lt. solver%max_qual)) then
!
                  n_qual_aop_in_sp  = n_qual_aop_in_sp + 1
                  n_qual_aop        = n_qual_aop + 1
!
               endif
!
            enddo
!
            if (n_qual_aop_in_sp .ne. 0) then
!
               n_qual_sp = n_qual_sp + 1
!
               call mem%alloc_int(sorted_qual_aop_in_sp_indices, n_qual_aop_in_sp, 1)
               call mem%alloc(sorted_qual_aop_in_sp, n_qual_aop_in_sp, 1)
!
               call get_n_highest(n_qual_aop_in_sp, last_sig_aop - first_sig_aop + 1, &
                                 D_xy(first_sig_aop:last_sig_aop, 1), sorted_qual_aop_in_sp, &
                                 sorted_qual_aop_in_sp_indices)
!
               n_previous_qual_aop = (n_qual_aop - n_qual_aop_in_sp)
!
               do aop = 1, n_qual_aop_in_sp
!
                  qual_aop(aop + n_previous_qual_aop, 1) = sig_aop_to_aos(sorted_qual_aop_in_sp_indices(aop, 1) &
                                                               + first_sig_aop - 1, 1)
!
                  qual_aop(aop + n_previous_qual_aop, 2) = sig_aop_to_aos(sorted_qual_aop_in_sp_indices(aop, 1) &
                                                               + first_sig_aop - 1, 2)
!
                  qual_aop(aop + n_previous_qual_aop, 3) = sorted_qual_aop_in_sp_indices(aop, 1) &
                                                               + first_sig_aop - 1
!
               enddo
!
               first_x = sig_aop_to_aos(first_sig_aop, 1) ! alpha
               first_y = sig_aop_to_aos(first_sig_aop, 2) ! beta
!
               qual_sp(n_qual_sp, 1) = system%basis2shell(first_x)
               qual_sp(n_qual_sp, 2) = system%basis2shell(first_y)
               qual_sp(n_qual_sp, 3) = n_qual_aop_in_sp
!
               call mem%dealloc_int(sorted_qual_aop_in_sp_indices, n_qual_aop_in_sp, 1)
               call mem%dealloc(sorted_qual_aop_in_sp, n_qual_aop_in_sp, 1)
!
            endif
!
            if (n_qual_aop == solver%max_qual) then
!
               exit
!
            endif
!
         enddo
!
         call mem%dealloc_int(sorted_max_sig_sp, n_sig_sp, 1)
!
!        Cut out the qualified parts of the aop and sp lists
!
         call mem%alloc_int(qual_aop_copy, n_qual_aop, 3)
         call mem%alloc_int(qual_sp_copy, n_qual_sp, 3)
!
         qual_aop_copy(:, :) = qual_aop(1 : n_qual_aop, :)
         qual_sp_copy(:, :)  = qual_sp(1 : n_qual_sp, :)
!
         call mem%dealloc_int(qual_aop, solver%max_qual, 3)
         call mem%dealloc_int(qual_sp, solver%n_sp, 3)
!
         call mem%alloc_int(qual_aop, n_qual_aop, 3)
         call mem%alloc_int(qual_sp, n_qual_sp, 3)
!
         qual_aop    = qual_aop_copy
         qual_sp     = qual_sp_copy
!
         call mem%dealloc_int(qual_aop_copy, n_qual_aop, 3)
         call mem%dealloc_int(qual_sp_copy, n_qual_sp, 3)
!
!        Prepare to construct g_wxyz in parallelized loop
!
         call mem%alloc_int(n_qual_aop_in_prev_sps, n_qual_sp, 1)
         n_qual_aop_in_prev_sps = 0
!
         do CD_sp = 1, n_qual_sp - 1
!
             n_qual_aop_in_prev_sps(CD_sp + 1, 1) = n_qual_aop_in_prev_sps(CD_sp, 1) + qual_sp(CD_sp, 3)
!
         enddo
!
!        Construct g_wxyz
!
         call mem%alloc(g_wxyz, n_sig_aop, n_qual_aop)
!
!$omp parallel do &
!$omp private(AB_sp, CD_sp, A, B, A_interval, B_interval, C, D, C_interval, D_interval, &
!$omp  aop, w, x, y, z, wx, yz, wx_packed, g_AB_CD, n_qual_aop_in_sp) &
!$omp shared(g_wxyz, n_qual_aop_in_prev_sps, qual_aop) &
!$omp schedule(guided)
         do CD_sp = 1, n_qual_sp
!
            C                = qual_sp(CD_sp, 1)
            D                = qual_sp(CD_sp, 2)
            n_qual_aop_in_sp = qual_sp(CD_sp, 3)
!
            C_interval = system%shell_limits(C)
            D_interval = system%shell_limits(D)
!
!           Calculate the ({wx} | J) integrals,
!           where {wx} is the screened list of integrals
!
            do AB_sp = 1, n_sig_sp
!
               A = sig_sp_to_shells(AB_sp, 1)
               B = sig_sp_to_shells(AB_sp, 2)
!
               A_interval = system%shell_limits(A)
               B_interval = system%shell_limits(B)
!
               call mem%alloc(g_AB_CD, &
                              (A_interval%size)*(B_interval%size), &
                              (C_interval%size)*(D_interval%size))
!
               call system%ao_integrals%get_ao_g_wxyz(g_AB_CD, A, B, C, D)
!
                 do aop = 1, n_qual_aop_in_sp
!
                    y = qual_aop(aop + n_qual_aop_in_prev_sps(CD_sp, 1), 1)
                    z = qual_aop(aop + n_qual_aop_in_prev_sps(CD_sp, 1), 2)
!
                    yz = C_interval%size*(z - D_interval%first) + y - C_interval%first + 1
!
                    if (A == B) then
!
                       do w = 1, A_interval%size
                          do x = w, B_interval%size
!
                             wx_packed = (max(w,x)*(max(w,x)-3)/2) + w + x
                             wx = A_interval%size*(x-1) + w
!
                             g_wxyz(sig_sp_to_first_sig_aop(AB_sp, 1) + wx_packed - 1, aop + n_qual_aop_in_prev_sps(CD_sp, 1)) &
                                   = g_AB_CD(wx, yz)
!
                          enddo
                       enddo
!
                    else
!
                       do w = 1, A_interval%size
                          do x = 1, B_interval%size
!
                             wx = A_interval%size*(x-1) + w
!
                             g_wxyz(sig_sp_to_first_sig_aop(AB_sp, 1) + wx - 1, aop + n_qual_aop_in_prev_sps(CD_sp, 1)) &
                                   = g_AB_CD(wx, yz)
!
                          enddo
                       enddo
!
                    endif
!
                 enddo
!
                  call mem%dealloc(g_AB_CD, &
                                    (A_interval%size)*(B_interval%size), &
                                    (C_interval%size)*(D_interval%size))
!
            enddo
!
         enddo ! cd_sp
!$omp end parallel do
!
         call mem%dealloc_int(n_qual_aop_in_prev_sps, n_qual_sp, 1)
!
!        Subtract old cholesky vectors
!
         call cpu_time(s_construct_time)
!
         if (solver%n_cholesky .ne. 0) then
!
            do I = 1, cholesky_array%n_nodes
!
               call cholesky_array%get_element(cholesky, I) ! Let cholesky point to node # I
!
               n_cholesky_in_node = cholesky_array%get_n_columns_element(I)
!
               call mem%alloc(cholesky_tmp, n_cholesky_in_node, n_qual_aop)
!
               do J = 1, n_qual_aop
!
                  cholesky_tmp(:, J) = cholesky(qual_aop(J, 3), :)
!
               enddo
!
               call dgemm('N', 'N',            &
                           n_sig_aop,          &
                           n_qual_aop,         &
                           n_cholesky_in_node, &
                           -one,               &
                           cholesky,           &
                           n_sig_aop,          &
                           cholesky_tmp,       &
                           n_cholesky_in_node, &
                           one,                &
                           g_wxyz,             &
                           n_sig_aop)
!
               call mem%dealloc(cholesky_tmp, n_cholesky_in_node, n_qual_aop)
!
            enddo
!
            call mem%alloc_int(cholesky_basis, solver%n_cholesky + n_qual_aop, 3)
            cholesky_basis(1 : solver%n_cholesky, :) = cholesky_basis_new(:, :)
            call mem%dealloc_int(cholesky_basis_new, solver%n_cholesky, 3)
!
         else
!
            call mem%alloc_int(cholesky_basis, n_qual_aop, 3)
            cholesky_basis = 0
!
         endif
!
         call mem%alloc(approx_diagonal_accumulative, n_sig_aop, 1)
         approx_diagonal_accumulative = zero
!
         call cholesky_array%push_back(n_sig_aop, n_qual_aop) ! Make new element at the tail
         call cholesky_array%get_element(cholesky_new, cholesky_array%n_nodes) ! Point cholesky_new to this last element
         cholesky_new = zero
!
         current_qual = 0
!
         construct_more_choleskys = .true.
!
         call mem%alloc_int(qual_max, n_qual_aop, 1)
!
         do while ((current_qual .lt. n_qual_aop) .and. construct_more_choleskys)
!
            current_qual = current_qual + 1
!
            D_max = zero
!
            do qual = 1, n_qual_aop
!
               xy = qual_aop(qual, 3)
!
               if (D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1) .gt. D_max) then
!
                  qual_max(current_qual, 1) = qual
                  D_max    = D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1)
                  xy_max = xy
!
               endif
!
            enddo
!
            if ((D_max*screening_vector(xy_max, 1) .gt. solver%threshold) .and. &
               (D_max .gt. solver%threshold)  .and. &
               (D_max .ge. solver%span*D_max_full)) then
!
               cholesky_basis(solver%n_cholesky + current_qual, 1) = qual_aop(qual_max(current_qual, 1), 1)
               cholesky_basis(solver%n_cholesky + current_qual, 2) = qual_aop(qual_max(current_qual, 1), 2)
!
               A = system%basis2shell(qual_aop(qual_max(current_qual, 1), 1))
               B = system%basis2shell(qual_aop(qual_max(current_qual, 1), 2))
!
               cholesky_basis(solver%n_cholesky + current_qual, 3) = get_sp_from_shells(A, B, solver%n_s)
!
               cholesky_new(: , current_qual) = g_wxyz(:, qual_max(current_qual, 1))
!
               if (current_qual .gt. 1) then
!
                  call mem%alloc(cholesky_tmp, 1, current_qual - 1)
!
                  cholesky_tmp(1, :) = cholesky_new(qual_aop(qual_max(current_qual, 1), 3), 1 : current_qual - 1)
!
                  call dgemm('N', 'T',                      &
                           n_sig_aop,                       &
                           1,                               &
                           current_qual - 1,                &
                           -one,                            &
                           cholesky_new,                    &
                           n_sig_aop,                       &
                           cholesky_tmp,                    &
                           1,                               &
                           one,                             &
                           cholesky_new(1, current_qual),   &
                           n_sig_aop)
!
                  call mem%dealloc(cholesky_tmp, 1, current_qual - 1)
!
               endif
!
               call dscal(n_sig_aop, one/sqrt(D_max), cholesky_new(1, current_qual), 1)
!
               do xy = 1, n_sig_aop
!
                  if (D_xy(xy, 1) == zero) then
!
                     cholesky_new(xy, current_qual) = zero
!
                  endif
!
               enddo
!
               do xy = 1, n_sig_aop
!
                  approx_diagonal_accumulative(xy, 1) = approx_diagonal_accumulative(xy, 1) + cholesky_new(xy, current_qual)**2
!
               enddo
!
               D_xy(qual_aop(qual_max(current_qual, 1), 3), 1) = zero
               approx_diagonal_accumulative(qual_aop(qual_max(current_qual, 1), 3), 1)  = zero
!
               do xy = 1, n_sig_aop
!
                  if (D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1) .lt. zero) then
!
                     if (abs(D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1)) .gt. 1.0d-10) then
                        if (write_warning) then
                           write(output%unit, '(a)') 'Warning: Found significant negative diagonal! '
                           write_warning = .false.
                        endif
                        sig_neg = sig_neg + 1
                     endif
!
                     D_xy(xy, 1) = zero
                     approx_diagonal_accumulative(xy, 1) = zero
!
                  elseif ((D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1))*screening_vector(xy,1) .lt. solver%threshold .or. &
                     (D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1)).lt. solver%threshold) then
!
                     D_xy(xy, 1) = zero
                     approx_diagonal_accumulative(xy, 1) = zero
!
                  endif
!
               enddo
!
            else
!
               construct_more_choleskys = .false.
               current_qual = current_qual - 1
!
            endif
!
         enddo ! End of decomposition
!
         do xy = 1, n_sig_aop
!
            D_xy(xy, 1) = D_xy(xy, 1) - approx_diagonal_accumulative(xy, 1)
!
         enddo
!
         call mem%dealloc(approx_diagonal_accumulative, n_sig_aop, 1)
!
         call mem%dealloc_int(qual_max, n_qual_aop, 1)
!
         call cpu_time(e_construct_time)
         full_construct_time = full_construct_time + e_construct_time - s_construct_time
!
         call mem%dealloc(g_wxyz, n_sig_aop, n_qual_aop)
!
         n_new_cholesky = current_qual
!
!        Find new significant diagonals
!
         n_new_sig_sp = 0
         allocate(new_sig_sp(n_sig_sp))
!
         new_sig_sp = .false.
!
         sig_sp_counter = 0
!
         do sp = 1, solver%n_sp
!
            if (sig_sp(sp)) then
!
               sig_sp_counter = sig_sp_counter + 1
!
               first = sig_sp_to_first_sig_aop(sig_sp_counter, 1)
               last  = sig_sp_to_first_sig_aop(sig_sp_counter + 1, 1) - 1
!
               new_sig_sp(sig_sp_counter) = (is_significant(D_xy(first:last, 1), &
                                                last - first + 1, solver%threshold, &
                                                screening_vector(first:last, 1) ) .and. &
                                             is_significant(D_xy(first:last, 1), &
                                                last - first + 1, solver%threshold ))
!
               sig_sp(sp) = new_sig_sp(sig_sp_counter)
!
               if (new_sig_sp(sig_sp_counter)) then
!
                  n_new_sig_sp = n_new_sig_sp + 1
!
               endif
!
            endif
!
         enddo
!
         call cpu_time(s_reduce_time)
!
         if (n_new_sig_sp .gt. 0) then
!
!           Update index lists: sps -> aops, aops -> aos, and sps -> full sps
!
            call mem%alloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
            new_sig_sp_to_first_sig_aop = 0
!
            call mem%alloc_int(sig_sp_to_previous_sig_sp, n_sig_sp + 1, 1) ! 1 2 3 4 ... n_sig_sp, n_sig_sp + 1
            sig_sp_to_previous_sig_sp(n_sig_sp + 1, 1) = n_sig_sp + 1
!
            current_new_sig_sp    = 1
            n_new_sig_aop = 0
            first_sig_aop = 1
!
            do sp = 1, n_sig_sp
!
               sig_sp_to_previous_sig_sp(sp, 1) = sp
!
               if (new_sig_sp(sp)) then
!
                  A = sig_sp_to_shells(sp, 1)
                  B = sig_sp_to_shells(sp, 2)
!
                  A_interval = system%shell_limits(A)
                  B_interval = system%shell_limits(B)
!
                  new_sig_sp_to_first_sig_aop(current_new_sig_sp, 1) = first_sig_aop
!
                  first_sig_aop = first_sig_aop + get_size_sp(A_interval, B_interval)
                  n_new_sig_aop = first_sig_aop - 1
!
                  current_new_sig_sp    = current_new_sig_sp + 1

               endif
!
            enddo
!
            new_sig_sp_to_first_sig_aop(current_new_sig_sp, 1) = n_new_sig_aop + 1
!
            call mem%alloc_int(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            call reduce_array_int(sig_aop_to_aos,       &
                               new_sig_aop_to_aos,      &
                               sig_sp_to_first_sig_aop, &
                               new_sig_sp,              &
                               n_sig_sp,                &
                               n_sig_aop,               &
                               n_new_sig_aop,           &
                               2)
!
            call mem%alloc_int(new_sig_sp_to_shells, n_new_sig_sp, 2)
            new_sig_sp_to_shells = 0
!
            call reduce_array_int(sig_sp_to_shells,        &
                                  new_sig_sp_to_shells,    &
                                  sig_sp_to_previous_sig_sp,  &
                                  new_sig_sp,              &
                                  n_sig_sp,                &
                                  n_sig_sp,                &
                                  n_new_sig_sp,            &
                                  2)
!
            call mem%dealloc_int(sig_sp_to_previous_sig_sp, n_sig_sp + 1, 1)
            call mem%dealloc_int(sig_sp_to_shells, n_sig_sp, 2)
            call mem%alloc_int(sig_sp_to_shells, n_new_sig_sp, 2)
!
            sig_sp_to_shells = new_sig_sp_to_shells
            call mem%dealloc_int(new_sig_sp_to_shells, n_new_sig_sp, 2)
!
            call mem%alloc(D_xy_new, n_new_sig_aop, 1)
!
           call reduce_vector(D_xy,                     &
                             D_xy_new,                  &
                             sig_sp_to_first_sig_aop,   &
                             new_sig_sp,                &
                             n_sig_sp,                  &
                             n_sig_aop,                 &
                             n_new_sig_aop)
!
            call mem%dealloc(D_xy, n_sig_aop, 1)
            call mem%alloc(D_xy, n_new_sig_aop, 1)
!
            call dcopy(n_new_sig_aop, D_xy_new, 1, D_xy, 1)
!
            call mem%dealloc(D_xy_new, n_new_sig_aop, 1)
!
            call mem%alloc(screening_vector_new, n_new_sig_aop, 1)
!
           call reduce_vector(screening_vector,         &
                             screening_vector_new,      &
                             sig_sp_to_first_sig_aop,   &
                             new_sig_sp,                &
                             n_sig_sp,                  &
                             n_sig_aop,                 &
                             n_new_sig_aop)
!
            call mem%dealloc(screening_vector, n_sig_aop, 1)
            call mem%alloc(screening_vector, n_new_sig_aop, 1)
!
            call dcopy(n_new_sig_aop, screening_vector_new, 1, screening_vector, 1)
!
            call mem%dealloc(screening_vector_new, n_new_sig_aop, 1)
!
!           Remove the unused columns of cholesky new
!
            call cholesky_array%keep_columns(cholesky_array%n_nodes, 1, n_new_cholesky)
!
            call cholesky_array%reduce(sig_sp_to_first_sig_aop,  &
                                       new_sig_sp,               &
                                       n_sig_sp,                 &
                                       n_new_sig_aop)
!
            cholesky_new => null()
!
!
            call mem%alloc_int(cholesky_basis_new, solver%n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : solver%n_cholesky + n_new_cholesky, :)
            call mem%dealloc_int(cholesky_basis, solver%n_cholesky + n_qual_aop, 3)
!
   !        Deallocate old lists & reallocate + copy over new lists
!
            deallocate(new_sig_sp)
!
            call mem%dealloc_int(sig_sp_to_first_sig_aop, n_sig_sp + 1, 1)
            call mem%alloc_int(sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
            sig_sp_to_first_sig_aop = new_sig_sp_to_first_sig_aop
            call mem%dealloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
!
            call mem%dealloc_int(sig_aop_to_aos, n_sig_aop, 2)
            call mem%alloc_int(sig_aop_to_aos, n_new_sig_aop, 2)
            sig_aop_to_aos = new_sig_aop_to_aos
            call mem%dealloc_int(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            n_sig_sp = n_new_sig_sp
            n_sig_aop = n_new_sig_aop
!
            solver%n_cholesky = solver%n_cholesky + n_new_cholesky
!
            call mem%dealloc_int(qual_aop, n_qual_aop, 3)
            call mem%dealloc_int(qual_sp, n_qual_sp, 3)
!
            write(output%unit, '(i4, 8x, i9, 1x, a1, i9, 6x, e12.5, 4x, i4, 8x, i7, 8x, i13)') &
            solver%iteration, n_sig_aop,'/', n_sig_sp, D_max_full , n_qual_aop, solver%n_cholesky, solver%n_cholesky*n_sig_aop
            flush(output%unit)
!
         else
!
            call cholesky_array%finalize()
!
            cholesky_new   => null()
            cholesky       => null()
!
            call mem%dealloc(D_xy, n_sig_aop, 1)
!
            call mem%alloc_int(cholesky_basis_new, solver%n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : solver%n_cholesky + n_new_cholesky, :)
            call mem%dealloc_int(cholesky_basis, solver%n_cholesky + n_qual_aop, 3)
!
            solver%n_cholesky = solver%n_cholesky + n_new_cholesky
!
            done = .true.
!
            write(output%unit, '(i4, 8x, i9, 1x, a1, i9, 6x, e12.5, 4x, i4, 8x, i7, 8x, i13)') &
            solver%iteration, 0,'/',0, D_max_full, n_qual_aop, solver%n_cholesky, 0
            flush(output%unit)
!
         endif
!
         call cpu_time(e_reduce_time)
         full_reduce_time = full_reduce_time + e_reduce_time - s_reduce_time
!
      enddo ! while not done
!
      write(output%unit, '(a)')&
       '----------------------------------------------------------------------------------------------------'
      flush(output%unit)
!
!     Timings
!
      call cpu_time(e_select_basis_time)
!
      write(output%unit,'(/a42, i7)')'Number of signigicant negative diagonals: ', sig_neg
!
!     Prepare info on basis
!
!     Construct a list of all shell pairs (and shells) that contain elements of the basis
!     and how many elements of the basis they contain
!
      call mem%alloc_int(basis_shell_info_full, solver%n_sp, 4) ! A, B, AB, n_basis_aops_in_sp
      basis_shell_info_full = 0
!
      n_sp_in_basis = 0
!
      do i = 1, solver%n_cholesky
!
         A = system%basis2shell(cholesky_basis_new(i, 1))
         B = system%basis2shell(cholesky_basis_new(i, 2))
!
         AB = get_sp_from_shells(A, B, solver%n_s)
!
         found = .false.
!
         do sp_in_basis = 1, n_sp_in_basis
!
            if (AB == basis_shell_info_full(sp_in_basis, 3)) then
               found = .true.
               basis_shell_info_full(sp_in_basis, 4) = basis_shell_info_full(sp_in_basis, 4) + 1
               exit
            endif
!
         enddo
!
         if(.not. found) then
!
            n_sp_in_basis = n_sp_in_basis + 1
!
            basis_shell_info_full(n_sp_in_basis, 1) = A
            basis_shell_info_full(n_sp_in_basis, 2) = B
            basis_shell_info_full(n_sp_in_basis, 3) = AB
            basis_shell_info_full(n_sp_in_basis, 4) = 1
!
         endif
!
      enddo
!
      call mem%alloc_int(basis_shell_info, n_sp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_sp_in_basis, :)
      call mem%dealloc_int(basis_shell_info_full, solver%n_sp, 4)
!
!     Write basis_shell_data file containing
!
!        1. number shell pairs in basis
!        2. basis_shell_info
!        3. cholesky_basis
!
      call disk%open_file(solver%basis_shell_data, 'write', 'rewind')
!
      write(solver%basis_shell_data%unit) n_sp_in_basis
      write(solver%basis_shell_data%unit) basis_shell_info
      write(solver%basis_shell_data%unit) cholesky_basis_new
!
      call disk%close_file(solver%basis_shell_data)
!
      call mem%dealloc_int(basis_shell_info, n_sp_in_basis, 4)
      call mem%dealloc_int(cholesky_basis_new, solver%n_cholesky, 3)
!
   end subroutine determine_auxilliary_cholesky_basis_eri_cd_solver
!
!
   subroutine construct_overlap_cholesky_vecs_eri_cd_solver(solver, system)
!!
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system) :: system
!
!     Local variables
!
!     Integers
!
      integer(i15) :: n_sp_in_basis, sp_in_basis
      integer(i15) :: n_vectors
      integer(i15) :: current_aop_in_sp
      integer(i15) :: A, B, C, D, AB, AB_sp, CD_sp
      integer(i15) :: I, J, K, L, KL
      integer(i15) :: w, x, y, z, wx, yz
!
!     Integer allocatable arrays
!
      integer(i15), dimension(:,:), allocatable :: basis_shell_info        ! Info on shells containing elements of the basis
      integer(i15), dimension(:,:), allocatable :: basis_shell_info_full   ! Info on shells containing elements of the basis
      integer(i15), dimension(:,:), allocatable :: cholesky_basis          ! ao and ao pair indices of the elements of the cholesky basis
      integer(i15), dimension(:,:), allocatable :: cholesky_basis_updated  ! ao and ao pair indices of the elements of the cholesky basis
      integer(i15), dimension(:,:), allocatable :: basis_aops_in_CD_sp     ! basis ao pairs in shell pair CD
      integer(i15), dimension(:,:), allocatable :: basis_aops_in_AB_sp     ! basis ao pairs in shell pair AB
!
      integer(kind=4), dimension(:), allocatable :: keep_vectors
!
!     Intervals
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
!     Reals
!
      real(dp) :: s_decomp_time, e_decomp_time, s_build_basis_time, e_build_basis_time
!
!     Real allocatable arrays
!
      real(dp), dimension(:,:), allocatable :: g_AB_CD
      real(dp), dimension(:,:), allocatable :: integrals_auxiliary, integrals_auxiliary_packed
      real(dp), dimension(:,:), allocatable :: cholesky_vecs
      real(dp), dimension(:,:), allocatable :: temp_cholesky
!
!     Logicals
!
      logical :: found
!
      call cpu_time(s_build_basis_time)
!
!     Read basis_shell_data
!
      call disk%open_file(solver%basis_shell_data, 'read')
      rewind(solver%basis_shell_data%unit)
!
      read(solver%basis_shell_data%unit) n_sp_in_basis
!
      call mem%alloc_int(basis_shell_info, n_sp_in_basis, 4)
      call mem%alloc_int(cholesky_basis, solver%n_cholesky, 3)
!
      read(solver%basis_shell_data%unit) basis_shell_info
      read(solver%basis_shell_data%unit) cholesky_basis
!
      call disk%close_file(solver%basis_shell_data, 'delete')
!
!
!     Construct integrals (J | J')
!
      call mem%alloc(integrals_auxiliary_packed, solver%n_cholesky*(solver%n_cholesky+1)/2, 1)
!
!$omp parallel do &
!$omp private(AB_sp, CD_sp, A, B, A_interval, B_interval, C, D, C_interval, D_interval, &
!$omp w, x, y, z, wx, yz, g_AB_CD, I, J, K, L, KL,&
!$omp current_aop_in_sp, basis_aops_in_CD_sp, basis_aops_in_AB_sp) &
!$omp shared(integrals_auxiliary_packed, cholesky_basis, basis_shell_info) &
!$omp schedule(guided)
      do CD_sp = 1, n_sp_in_basis
!
         C = basis_shell_info(CD_sp, 1)
         D = basis_shell_info(CD_sp, 2)
!
         C_interval = system%shell_limits(C)
         D_interval = system%shell_limits(D)
!
         call mem%alloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
!        Determine which elements in the shell pair CD are elements of the basis
!
         current_aop_in_sp = 0
!
         do I = 1, solver%n_cholesky
            if (cholesky_basis(I,3) == basis_shell_info(CD_sp, 3)) then
!
               current_aop_in_sp = current_aop_in_sp + 1
!
               basis_aops_in_CD_sp(current_aop_in_sp, 1) = cholesky_basis(I,1) - C_interval%first + 1
               basis_aops_in_CD_sp(current_aop_in_sp, 2) = cholesky_basis(I,2) - D_interval%first + 1
               basis_aops_in_CD_sp(current_aop_in_sp, 3) = I
!
            endif
         enddo
!
         do AB_sp = 1, n_sp_in_basis
!
            A = basis_shell_info(AB_sp, 1)
            B = basis_shell_info(AB_sp, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            call mem%alloc_int(basis_aops_in_AB_sp, basis_shell_info(AB_sp, 4), 3)
!
!           Determine which elements in the shell pair AB are elements of the basis
!
            current_aop_in_sp = 0
!
            do I = 1, solver%n_cholesky
               if (cholesky_basis(I,3) == basis_shell_info(AB_sp, 3)) then
!
                  current_aop_in_sp = current_aop_in_sp + 1
!
                  basis_aops_in_AB_sp(current_aop_in_sp, 1) = cholesky_basis(I,1) - A_interval%first + 1
                  basis_aops_in_AB_sp(current_aop_in_sp, 2) = cholesky_basis(I,2) - B_interval%first + 1
                  basis_aops_in_AB_sp(current_aop_in_sp, 3) = I
!
               endif
            enddo
!
!           Construct integrals
!
            call mem%alloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
            call system%ao_integrals%get_ao_g_wxyz(g_AB_CD, A, B, C, D)
!
!           Only keep those that correspond to elements of the basis
!
            do I = 1, basis_shell_info(AB_sp, 4)
               do J = 1, basis_shell_info(CD_sp, 4)
!
                  y = basis_aops_in_CD_sp(J, 1)
                  z = basis_aops_in_CD_sp(J, 2)
                  yz = C_interval%size*(z-1)+y
!
                  w = basis_aops_in_AB_sp(I, 1)
                  x = basis_aops_in_AB_sp(I, 2)
                  wx = A_interval%size*(x-1)+w
!
                  K = basis_aops_in_AB_sp(I, 3)
                  L = basis_aops_in_CD_sp(J, 3)
                  KL = (max(K,L)*(max(K,L)-3)/2) + K + L
!
                  integrals_auxiliary_packed(KL, 1) = g_AB_CD(wx, yz)
!
               enddo
            enddo
!
            call mem%dealloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
            call mem%dealloc_int(basis_aops_in_AB_sp, basis_shell_info(AB_sp, 4), 3)
!
         enddo ! AB
!
         call mem%dealloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
      enddo ! CD
!$omp end parallel do
!
      call mem%dealloc_int(basis_shell_info, n_sp_in_basis, 4)
!
!     Square up integrals
!
      call mem%alloc(integrals_auxiliary, solver%n_cholesky, solver%n_cholesky)
      call squareup(integrals_auxiliary_packed, integrals_auxiliary, solver%n_cholesky)
      call mem%dealloc(integrals_auxiliary_packed, solver%n_cholesky*(solver%n_cholesky + 1)/2, 1)
!
      n_vectors = 0
      allocate(keep_vectors(solver%n_cholesky))
!
      call cpu_time(s_decomp_time)
!
      call mem%alloc(temp_cholesky, solver%n_cholesky, solver%n_cholesky)
!
      call full_cholesky_decomposition_system(integrals_auxiliary, temp_cholesky, solver%n_cholesky, n_vectors, &
                                                solver%threshold*1.0d-1, keep_vectors)
!
      call cpu_time(e_decomp_time)
!
      call mem%dealloc(integrals_auxiliary, solver%n_cholesky, solver%n_cholesky)
!
      call mem%alloc(cholesky_vecs, n_vectors, n_vectors)
!
      do  J = 1, n_vectors
         do I = 1, n_vectors
!
            cholesky_vecs(I, J) = temp_cholesky(I, J)
!
         enddo
      enddo
!
      call mem%dealloc(temp_cholesky, solver%n_cholesky, solver%n_cholesky)
!
      call mem%alloc_int(cholesky_basis_updated, n_vectors, 3)
!
      do I = 1, n_vectors
!
         cholesky_basis_updated(I, :) = cholesky_basis(keep_vectors(I), :)
!
      enddo
!
      call mem%dealloc_int(cholesky_basis, solver%n_cholesky, 3)
      deallocate(keep_vectors)
!
      solver%n_cholesky = n_vectors
!
      write(output%unit, '(t6, a34, i7)')'Final number of cholesky vectors: ', solver%n_cholesky
      flush(output%unit)
!
!     Update the basis_shell_info array which contains information of which shell pairs (and shells)
!     contain elements of the basis and how many elements of the basis they contain.
!
      call mem%alloc_int(basis_shell_info_full, solver%n_sp, 4) ! A, B, AB, n_basis_aops_in_sp
      basis_shell_info_full = 0
!
      n_sp_in_basis = 0
!
      do i = 1, solver%n_cholesky
!
         A = system%basis2shell(cholesky_basis_updated(i, 1))
         B = system%basis2shell(cholesky_basis_updated(i, 2))
!
         AB = get_sp_from_shells(A, B, solver%n_s)
!
         found = .false.
!
         do sp_in_basis = 1, n_sp_in_basis
!
              if (AB == basis_shell_info_full(sp_in_basis, 3)) then
                 found = .true.
                 basis_shell_info_full(sp_in_basis, 4) = basis_shell_info_full(sp_in_basis, 4) + 1
                 exit
              endif
!
         enddo
!
         if(.not. found) then
!
            n_sp_in_basis = n_sp_in_basis + 1
!
            basis_shell_info_full(n_sp_in_basis, 1) = A
            basis_shell_info_full(n_sp_in_basis, 2) = B
            basis_shell_info_full(n_sp_in_basis, 3) = AB
            basis_shell_info_full(n_sp_in_basis, 4) = 1
!
         endif
!
      enddo
!
      call mem%alloc_int(basis_shell_info, n_sp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_sp_in_basis, :)
      call mem%dealloc_int(basis_shell_info_full, solver%n_sp, 4)
!
!     Write basis_shell_data file containing
!
!        1. number shell pairs in basis
!        2. basis_shell_info
!        3. cholesky_basis
!
      call disk%open_file(solver%basis_shell_data, 'write')
!
      write(solver%basis_shell_data%unit) n_sp_in_basis
      write(solver%basis_shell_data%unit) basis_shell_info
      write(solver%basis_shell_data%unit) cholesky_basis_updated
!
      call disk%close_file(solver%basis_shell_data)
!
      call mem%dealloc_int(basis_shell_info, n_sp_in_basis, 4)
      call mem%dealloc_int(cholesky_basis_updated, n_vectors, 3)
!
!     Write cholesky_aux file containing
!
!        1. Cholesky vectors L_JK
!
      call disk%open_file(solver%cholesky_aux, 'write', 'rewind')
!
      write(solver%cholesky_aux%unit) cholesky_vecs
!
      call disk%close_file(solver%cholesky_aux)
!
      call mem%dealloc(cholesky_vecs, n_vectors, n_vectors)
!
      call cpu_time(e_build_basis_time)
!
   end subroutine construct_overlap_cholesky_vecs_eri_cd_solver
!
!
   subroutine construct_cholesky_vectors_eri_cd_solver(solver, system)
!!
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      type(molecular_system) :: system
!
!     Local variables
!
!     Integers
!
      integer(i15) :: A, B, AB_sp, C, D, CD_sp
      integer(i15) :: w, x, y, z, wx, yz, yz_packed
      integer(i15) :: L, J, I
      integer(i15) :: n_sig_sp, n_sig_aop
      integer(i15) :: n_sp_in_basis, last_sp_included, sp_counter
      integer(i15) :: current_aop_in_sp
      integer(i15) :: n_AB_included, n_AB_included_current
      integer(i15) :: rec_offset
      integer(i15) :: size_AB, size_AB_current
!
!     Integer allocatable arrays
!
      integer(i15), dimension(:,:), allocatable :: basis_shell_info     ! Info on shells containing elements of the basis
      integer(i15), dimension(:,:), allocatable :: AB_info              ! Info on offsets and shells for OMP-loop [offset, A, B]
      integer(i15), dimension(:,:), allocatable :: basis_aops_in_CD_sp  ! Basis ao pairs in shell pair CD
      integer(i15), dimension(:,:), allocatable :: basis_aops_in_AB_sp  ! Basis ao pairs in shell pair AB
      integer(i15), dimension(:,:), allocatable :: cholesky_basis       ! Info on cholesky basis
!
!     Reals
!
      real(dp) :: s_construct_time, e_construct_time, full_construct_time
      real(dp) :: s_build_vectors_time, e_build_vectors_time
!
!     Real allocatable arrays
!
      real(dp), dimension(:,:), allocatable :: g_J_yz, g_CD_AB
      real(dp), dimension(:,:), allocatable :: L_K_yz
      real(dp), dimension(:,:), allocatable :: aux_chol_inverse
!
!     Logicals
!
      logical :: done, found_size
!
!     Logical allocatable arrays
!
      logical, dimension(:), allocatable :: sig_sp
!
!     Intervals
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
      write(output%unit, '(/a)') '- Construct Cholesky vectors'
      flush(output%unit)
!
      call cpu_time(s_build_vectors_time)
!
!     Read diagonal info
!
      call disk%open_file(solver%diagonal_info_target, 'read')
!
      allocate(sig_sp(solver%n_sp))
!
      read(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
      read(solver%diagonal_info_target%unit) sig_sp
!
      call disk%close_file(solver%diagonal_info_target)
!
!     Read inverse of cholesky vectors of auxiliary overlap
!
      call disk%open_file(solver%cholesky_aux_inverse, 'read')
!
      call mem%alloc(aux_chol_inverse, solver%n_cholesky, solver%n_cholesky)
!
      read(solver%cholesky_aux_inverse%unit) aux_chol_inverse
!
      call disk%close_file(solver%cholesky_aux_inverse)
!
!     Prepare file for AO Cholesky vectors
!
      call solver%cholesky_ao_vectors%init('cholesky_ao', 'sequential', 'unformatted')
      call disk%open_file(solver%cholesky_ao_vectors, 'write', 'rewind')

      call solver%cholesky_ao_vectors_info%init('cholesky_ao_batching_info', 'sequential', 'formatted')
      call disk%open_file(solver%cholesky_ao_vectors_info, 'write', 'rewind')
!
      rec_offset = 0
!
!     Construct (K | yz) and do matrix multiplication
!     sum_K (K | J)^-1 (J | yz) in batches of yz
!
      done = .false.
!
      full_construct_time = zero
!
      do while (.not. done)
!
!        Determine size of batch
!
         sp_counter              = 0
         n_AB_included_current   = 0
         size_AB_current         = 0
!
         found_size = .false.
!
         do B = 1, solver%n_s
            do A = B, solver%n_s
!
               sp_counter = sp_counter + 1
!
               A_interval = system%shell_limits(A)
               B_interval = system%shell_limits(B)
!
               if (sig_sp(sp_counter)) then
!
                  size_AB_current = size_AB_current + get_size_sp(A_interval, B_interval)
                  n_AB_included_current = n_AB_included_current + 1
!
                  if ((2*size_AB_current*(solver%n_cholesky)*dp + (solver%n_cholesky**2)*dp)*1.1d0 .ge. mem%available) then ! 10 percent buffer
!
                     if (.not. found_size) then
!
                        size_AB = size_AB_current - get_size_sp(A_interval, B_interval)
!
                        last_sp_included  = sp_counter - 1
                        n_AB_included     = n_AB_included_current - 1
!
                        found_size = .true.
!
                     endif
!
                  endif
!
               endif
!
            enddo
         enddo
!
         if (.not. found_size) then
!
            size_AB          = size_AB_current
            last_sp_included = sp_counter
            n_AB_included    = n_AB_included_current
!
         endif
!
         write(solver%cholesky_ao_vectors_info%unit, *) size_AB
!
         call disk%open_file(solver%basis_shell_data, 'read')
!
         read(solver%basis_shell_data%unit) n_sp_in_basis
!
         call mem%alloc_int(basis_shell_info, n_sp_in_basis, 4)
         call mem%alloc_int(cholesky_basis, solver%n_cholesky, 3)
!
         read(solver%basis_shell_data%unit) basis_shell_info
         read(solver%basis_shell_data%unit) cholesky_basis
!
         call disk%close_file(solver%basis_shell_data)
!
!        Construct g_J_yz = (J | yz)
!
         call mem%alloc(g_J_yz, solver%n_cholesky, size_AB)
         call mem%alloc_int(AB_info, n_AB_included, 3) ! [offset, A, B]
         AB_info = zero
!
         sp_counter = 0
!
         do B = 1, solver%n_s
            do A = B, solver%n_s
!
               AB_sp = get_sp_from_shells(A, B, solver%n_s)
!
               if (sig_sp(AB_sp) .and. AB_sp .le. last_sp_included) then
!
                  sp_counter = sp_counter + 1
                  sig_sp(AB_sp) = .false.
!
                  A_interval = system%shell_limits(A)
                  B_interval = system%shell_limits(B)
!
                  if (sp_counter .lt. n_AB_included) AB_info(sp_counter + 1, 1) = AB_info(sp_counter, 1) &
                                                + get_size_sp(A_interval, B_interval)
!
                  AB_info(sp_counter, 2) = A
                  AB_info(sp_counter, 3) = B
!
               endif
!
            enddo
         enddo
!
!$omp parallel do &
!$omp private(AB_sp, CD_sp, I, A, B, A_interval, &
!$omp B_interval, C, D, C_interval, D_interval, &
!$omp basis_aops_in_CD_sp, current_aop_in_sp, g_CD_AB, &
!$omp w, x, y, z, wx, yz, yz_packed, L, J) &
!$omp shared(g_J_yz, AB_info, basis_shell_info, cholesky_basis) &
!$omp schedule(guided)
         do AB_sp = 1, n_AB_included
!
            A = AB_info(AB_sp, 2)
            B = AB_info(AB_sp, 3)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            do CD_sp = 1, n_sp_in_basis
!
               C = basis_shell_info(CD_sp, 1)
               D = basis_shell_info(CD_sp, 2)
!
               C_interval = system%shell_limits(C)
               D_interval = system%shell_limits(D)
!
               call mem%alloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
               current_aop_in_sp = 0
!
               do I = 1, solver%n_cholesky
                  if (cholesky_basis(I,3) == basis_shell_info(CD_sp, 3)) then
!
                     current_aop_in_sp = current_aop_in_sp + 1
!
                     basis_aops_in_CD_sp(current_aop_in_sp, 1) = cholesky_basis(I,1) - C_interval%first + 1
                     basis_aops_in_CD_sp(current_aop_in_sp, 2) = cholesky_basis(I,2) - D_interval%first + 1
                     basis_aops_in_CD_sp(current_aop_in_sp, 3) = I
!
                  endif
               enddo
!
              call mem%alloc(g_CD_AB, &
                       (C_interval%size)*(D_interval%size), &
                       (A_interval%size)*(B_interval%size))
!
              g_CD_AB = zero
!
              call system%ao_integrals%get_ao_g_wxyz(g_CD_AB, C, D, A, B)
!
               if (A == B) then
!
                  do y = 1, A_interval%size
                     do z = y, B_interval%size
!
                        yz_packed = (max(y,z)*(max(y,z)-3)/2) + y + z
                        yz = A_interval%size*(z-1) + y
!
                        do J = 1, basis_shell_info(CD_sp, 4)
                           w = basis_aops_in_CD_sp(J, 1)
                           x = basis_aops_in_CD_sp(J, 2)
                           L = basis_aops_in_CD_sp(J, 3)
                           wx = C_interval%size*(x-1)+w

                           g_J_yz(L, yz_packed + AB_info(AB_sp, 1)) = g_CD_AB(wx, yz)
!
                           enddo
                        enddo
                     enddo
!
                  else
!
                     do y = 1, A_interval%size
                        do z = 1, B_interval%size
                           do J = 1, basis_shell_info(CD_sp, 4)
                              w = basis_aops_in_CD_sp(J, 1)
                              x = basis_aops_in_CD_sp(J, 2)
                              L = basis_aops_in_CD_sp(J, 3)
!
                              wx = C_interval%size*(x-1) + w
!
                              yz = A_interval%size*(z-1) + y
!
                              g_J_yz(L, yz + AB_info(AB_sp, 1)) = g_CD_AB(wx, yz)
!
                           enddo
                        enddo
                     enddo
!
                  endif
!
                  call mem%dealloc(g_CD_AB,                 &
                     (C_interval%size)*(D_interval%size),   &
                     (A_interval%size)*(B_interval%size))
!
               call mem%dealloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
            enddo ! CD
!
         enddo ! AB
!$omp end parallel do
!
         call mem%dealloc_int(AB_info, n_AB_included, 3)
!
         call cpu_time(s_construct_time)
!
!        L_K_yz = sum_J (K | J)^-1 (J | yz)
!
         call mem%alloc(L_K_yz, solver%n_cholesky, size_AB)
!
         call dgemm('N', 'N',             &
                       solver%n_cholesky, &
                       size_AB,           &
                       solver%n_cholesky, &
                       one,               &
                       aux_chol_inverse,  & !(K|J)^-1
                       solver%n_cholesky, &
                       g_J_yz,            &
                       solver%n_cholesky, &
                       zero,              &
                       L_K_yz,            &
                       solver%n_cholesky)
!
         call mem%dealloc(g_J_yz, solver%n_cholesky, size_AB)
!
         call cpu_time(e_construct_time)
         full_construct_time = full_construct_time + e_construct_time - s_construct_time
!
!        Write vectors to file
!
         do J = 1, solver%n_cholesky
!
            write(solver%cholesky_ao_vectors%unit) (L_K_yz(J, I), I = 1, size_AB)
!
         enddo
!
         rec_offset = rec_offset + size_AB
!
         call mem%dealloc(L_K_yz, solver%n_cholesky, size_AB)
!
         done = .true.
!
         do I = 1, solver%n_sp
            if (sig_sp(I)) then
!
               done = .false.
               exit
!
            endif
         enddo
!
      enddo ! done
!
      call disk%close_file(solver%cholesky_ao_vectors)
!
      write(solver%cholesky_ao_vectors_info%unit, *) 'DONE'
!
      call disk%close_file(solver%cholesky_ao_vectors_info)
!
      call mem%dealloc(aux_chol_inverse, solver%n_cholesky, solver%n_cholesky)
      call mem%dealloc_int(cholesky_basis, solver%n_cholesky, 3)
      call mem%dealloc_int(basis_shell_info, n_sp_in_basis, 4)
      deallocate(sig_sp)
!
   end subroutine construct_cholesky_vectors_eri_cd_solver
!
!
   subroutine invert_overlap_cholesky_vecs_eri_cd_solver(solver)
!!
!!    Invert cholesky vectors of auxiliary basis overlap
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, July 2018
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(file) :: cholesky_file, cholesky_inverse_file
!
      real(dp):: s_invert_time, e_invert_time
!
      real(dp), dimension(:,:), allocatable :: cholesky, cholesky_inverse
!
      write(output%unit, '(/a)') '- Inverting auxiliary basis'
      flush(output%unit)
      call cpu_time(s_invert_time)
!
!     Read Cholesky vectors of auxiliary basis overlap
!
      call disk%open_file(solver%cholesky_aux, 'read')
      rewind(solver%cholesky_aux%unit)
!
      call mem%alloc(cholesky, solver%n_cholesky, solver%n_cholesky)
!
      read(solver%cholesky_aux%unit) cholesky
!
      call disk%close_file(solver%cholesky_aux)
!
!     Invert cholesky vectors
!
      call mem%alloc(cholesky_inverse, solver%n_cholesky, solver%n_cholesky)
!
      call inv_lower_tri(cholesky_inverse, cholesky, solver%n_cholesky)
!
      call mem%dealloc(cholesky, solver%n_cholesky, solver%n_cholesky)
!
!     Write inverse Cholesky vectors of auxiliary basis overlap
!
!        1. n_cholesky: Number of elements of the basis
!        2. cholesky_inverse (n_cholesky, n_cholesky)
!
      call disk%open_file(solver%cholesky_aux_inverse, 'write', 'rewind')
!
      write(solver%cholesky_aux_inverse%unit) cholesky_inverse
!
      call disk%close_file(solver%cholesky_aux_inverse)
!
      call mem%dealloc(cholesky_inverse, solver%n_cholesky, solver%n_cholesky)
!
!     Timings
!
      call cpu_time(e_invert_time)
!
!
   end subroutine invert_overlap_cholesky_vecs_eri_cd_solver
!
!
   subroutine cholesky_vecs_diagonal_test_eri_cd_solver(solver)
!!
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(file) :: cholesky_ao_vectors
!
      real(dp), dimension(:,:), allocatable :: D_diff, L_K_yz
!
      real(dp) :: ddot, max_diff, min_diff
!
      integer(i15) :: aop, n_sig_sp, n_sig_aop, J, I, size_AB, AB_offset
!
      character(len=40) :: line
!
      write(output%unit, '(/a)') '- Test Cholesky vectors NOT WORKING!!'
      flush(output%unit)
!
!     Read diagonal information
!
      call disk%open_file(solver%diagonal_info_target, 'read')
!
      rewind(solver%diagonal_info_target)
!
      read(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
!
      call mem%alloc(D_diff, n_sig_aop, 1)
!
      read(solver%diagonal_info_target%unit)
      read(solver%diagonal_info_target%unit) D_diff
!
      call disk%close_file(solver%diagonal_info_target)
!
      call disk%open_file(solver%cholesky_ao_vectors, 'read')
      call disk%open_file(solver%cholesky_ao_vectors_info, 'read')
!
      rewind(solver%cholesky_ao_vectors%unit)
      rewind(solver%cholesky_ao_vectors_info%unit)
!
      read(solver%cholesky_ao_vectors_info%unit, *) line
!
      AB_offset = 0
!
      do while (trim(line) .ne. 'DONE')
!
!        Calculate difference between actual and approximate diagonal
!
         read(line, *) size_AB
!
         call mem%alloc(L_K_yz, solver%n_cholesky, size_AB)
!
         do J = 1, solver%n_cholesky
!
            read(solver%cholesky_ao_vectors%unit) (L_K_yz(J, I), I = 1, size_AB)
!
         enddo
!
         do I = 1, size_AB 
!
            D_diff(I + AB_offset, 1) = D_diff(I + AB_offset, 1) - ddot(solver%n_cholesky, L_K_yz(1, I), 1, L_K_yz(1, I), 1)
!
         enddo
!
         call mem%dealloc(L_K_yz, solver%n_cholesky, size_AB)
!
         AB_offset = AB_offset + size_AB
         read(solver%cholesky_ao_vectors_info%unit, *) line
!
      enddo
!
      call disk%close_file(solver%cholesky_ao_vectors)
      call disk%close_file(solver%cholesky_ao_vectors_info)
!
!     Calculate maximal difference and minimal difference
!
      max_diff = zero
!
      do aop = 1, n_sig_aop
         if (abs(D_diff(aop, 1)) .gt. max_diff) max_diff = abs(D_diff(aop, 1))
      enddo
!
      min_diff =1.0d5
!
      do aop = 1, n_sig_aop
         if (D_diff(aop, 1) .lt. min_diff) min_diff = D_diff(aop, 1)
      enddo
!
      call mem%dealloc(D_diff, n_sig_aop, 1)
!
      write(output%unit, '(/a71, e12.4)')'Maximal difference between approximate and actual diagonal:            ', max_diff
      write(output%unit, '(/a71, e12.4)')'Minimal element of difference between approximate and actual diagonal: ', min_diff
      flush(output%unit)
!
   end subroutine cholesky_vecs_diagonal_test_eri_cd_solver
!
!
   integer(i15) function get_size_sp(A_interval, B_interval)
!!
!!    Get size shell pair
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns size of diagonal for a given shell pair.
!!
      implicit none
!
      type(interval) :: A_interval
      type(interval) :: B_interval
!
      if (A_interval%first == B_interval%first) then
!
         get_size_sp = A_interval%size*(A_interval%size + 1)/2
!
      else
!
         get_size_sp = (A_interval%size)*(B_interval%size)
!
      endif
!
   end function get_size_sp
!
!
   integer(i15) function get_sp_from_shells(s1, s2, n_s)
!!
!!    Get shell pair from shells,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      integer(i15), intent(in) :: s1, s2, n_s
      integer(i15) :: A, B
!
      get_sp_from_shells = 1
!
      do B = 1, n_s
         do A = B, n_s
!
            if (s1 == A .and. s2 == B)  return
!
            get_sp_from_shells = get_sp_from_shells + 1
!
         enddo
      enddo

!
   end function get_sp_from_shells
!
!
 subroutine read_info_eri_cd_solver(solver)
!!
!!    Read information
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!     Read input if it is present:
!!
!!        cholesky
!!           threshold: 1.0d-8
!!           span: 1.0d-2
!!           qualified: 1000
!!           one center
!!           no vectors
!!        end cholesky
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      character(len=100) :: line
      character(len=100) :: current_basis
!
      integer(i15) :: i = 0, ioerror
!
      rewind(input%unit)
!
      read(input%unit,'(a)') line
      line = remove_preceding_blanks(line)
!
      do while ((trim(line) .ne. 'end cholesky') .and. (line(1:2) .ne. 'geometry'))
!
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)

         if (trim(line) == 'cholesky') then ! found cholesky section in input
!
            do while (trim(line) .ne. 'end cholesky')
!
               read(input%unit,'(a)') line
               line = remove_preceding_blanks(line)
!
               if (line(1:10) == 'threshold:') then
!
                  read(line(11:100), '(d16.5)') solver%threshold
!
               elseif (line(1:5) == 'span:') then
!
                  read(line(6:100), '(d16.5)') solver%span
!
               elseif (line(1:10) == 'qualified:') then
!
                  read(line(11:100), '(d16.5)') solver%max_qual
!
               elseif (trim(line) == 'one center') then
!
                  solver%one_center = .true.
!
               elseif (trim(line) == 'no vectors') then
!
                  solver%construct_vectors = .false.
!
               endif
!
            end do
!
            backspace(input%unit)
!
         endif
!
      enddo
      backspace(input%unit)
!
   end subroutine read_info_eri_cd_solver
!
!
   subroutine construct_mo_cholesky_vecs(solver, system, n_mo, orbital_coefficients)
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      call disk%open_file(solver%cholesky_ao_vectors, 'read')
      call disk%open_file(cholesky_mo_vectors_seq, 'write', 'rewind')
      call disk%open_file(solver%cholesky_mo_vectors, 'write', 'rewind')
      call disk%open_file(solver%cholesky_ao_vectors_info, 'read')
!
!     Read significant sp info 
!
!
      call disk%open_file(solver%diagonal_info_target, 'read')
!
      rewind(solver%diagonal_info_target)
!
      read(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
!
      allocate(sig_sp(solver%n_s))
!
      read(solver%diagonal_info_target%unit) sig_sp
!
      call disk%close_file(solver%diagonal_info_target)
!
     do J = 1, solver%n_cholesky
!
        rewind(solver%cholesky_ao_vectors%unit)
        rewind(solver%cholesky_ao_vectors_info%unit)
!
        read(solver%cholesky_ao_vectors_info%unit, *) line
!
        AB_offset = 0
!
         call mem%alloc(L_J_xy, 1, n_sig_aop)
         call mem%alloc(L_J_yx_full, solver%n_ao, solver%n_ao)
!
         L_J_xy = zero
         L_J_yx_full = zero
!
         do while (trim(line) .ne. 'DONE')
!
!          Calculate difference between actual and approximate diagonal
!
           read(line, *) size_AB     
!
!          Empty reads 
!
           do I = 1, J - 1
!
              read(solver%cholesky_ao_vectors%unit)
!
           enddo
!
           read(solver%cholesky_ao_vectors%unit) (L_J_yz(1, AB_offset + I), I = 1, size_AB)
!
!          Empty reads 
!
            do I = J + 1, solver%n_cholesky 
!
              read(solver%cholesky_ao_vectors%unit)
!
           enddo
!
           AB_offset = AB_offset + size_AB
!
           read(solver%cholesky_ao_vectors_info%unit, *) line
!
        enddo
!
         sp = 0
!
         do B = 1, solver%n_s
!
            B_interval = system%shell_limits(A)
!
            do A = B, solver%n_s
!
               A_interval = system%shell_limits(B)
!
               sp = sp + 1
!
               if (sig_sp(sp)) then
!
                  if (A .ne. B) then 
!
                     do x = 1, A_interval%size
                        do y = 1, B_interval%size
!  
                            xy = A_interval%size*(y-1) + x
                            L_J_xy_full(1, xy + AB_offset) = L_J_xy(1, xy + AB_offset)
!  
                        enddo
                     enddo
!
                  else
!
                     do x = 1, A_interval%size
                        do y = 1, A_interval%size
!  
                            xy = A_interval%size*(y-1) + x
                            xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
                            L_J_xy_full(1, xy + AB_offset) = L_J_xy(1, xy + AB_offset)
!  
                        enddo
                     enddo
!
               else
!
               AB_offset_full = AB_offset_full + (A_interval%size)*(B_interval%size)
               AB_offset = AB_offset + get_size_sp(A_interval, B_interval)
!
            enddo
!
            call mem%dealloc(L_J_xy, 1, n_sig_aop)
            call mem%dealloc(L_J_yx_full, solver%n_ao, solver%n_ao)
!
         enddo
!
!       Transform the AO vectors to form the Cholesky MO vectors
!
        call mem%alloc(X, n_ao, n_mo)
!
        call dgemm('N','N',               &
                    solver%n_ao,          &
                    n_mo,                 &
                    solver%n_ao,          &
                    one,                  &
                    L_J_xy_full,          &
                    n_ao,                 &
                    orbital_coefficients, &
                    solver%n_ao,          &
                    zero,                 &
                    X,                    &
                    solver%n_ao)
!
         call mem%dealloc(L_J_pq, n_mo, n_mo)
!
        call dgemm('T','N',               &
                    n_mo,                 &
                    n_mo,                 &
                    solver%n_ao,          &
                    one,                  &
                    orbital_coefficients, &
                    solver%n_ao,          &
                    X,                    &
                    solver%n_ao,          &
                    zero,                 &
                    L_J_pq,               &
                    n_mo)
!
         call mem%dealloc(X, solver%n_ao, n_mo)
!
        write(cholesky_mo_vectors_seq%unit) ((L_J_pq(p,q), q = 1, p), q = 1, n_mo)
        call mem%dealloc(L_J_pq, n_mo, n_mo)
!
      enddo
!
!     Read L_pq_J in batches over q
!
      required = ((n_mo)**2)*(solver%n_cholesky)
!
!     Initialize batching variable
!
      call batch_q%init(n_mo)
      call mem%num_batch(batch_q, required)
!
!     Loop over the q-batches
!
      do current_q_batch = 1, batch_q%num_batches
!
!        Determine limits of the q-batch
!
         call batch_q%determine_limits(current_q_batch)
!
         call mem%alloc(L_pq_J,                                            &
            (((batch_q%length + 1)*(batch_q%length)/2) +                   &
            (n_mo - batch_q%length - batch_q%first + 1)*(batch_q%length)), &
            solver%n_cholesky)
!
         if (batch_q%first .ne. 1) then
!
!           Calculate index of last element to throw away
!
            throw_away_index = index_packed(n_mo, batch_q%first - 1)
!
!           Throw away all elements from 1 to throw_away_index, then read from batch start
!
            do j = 1, solver%n_cholesky
!
              read(cholesky_mo_vectors_seq%unit) (throw_away, i = 1, throw_away_index), &
                  (L_pq_J(p,j), p = 1, (((batch_q%length + 1)*(batch_q%length)/2) + &
                                       (n_mo - batch_q%length - batch_q%first + 1)*(batch_q%length)))
!
            enddo
!
         else
!
!           Read from the start of each entry
!
            do j = 1, solver%n_cholesky
!
              read(cholesky_mo_vectors_seq%unit) (L_pq_J(p,j), p = 1, (((batch_q%length + 1)*(batch_q%length)/2) + &
                                    (n_mo - batch_q%length - batch_q%first + 1)*(batch_q%length)))
!
            enddo
!
         endif
!
         do p = 1, n_mo
            do q = batch_q%first, batch_q%last
!
               pq = (max(p, q)*(max(p, q)-3)/2) + p + q
               write(solver%cholesky_mo_vectors%unit, rec=pq) (L_pq_J(pq, J), J = 1, solver%n_cholesky)
!
            enddo
         enddo
!
         call mem%dealloc(L_pq_J,                                          &
            (((batch_q%length + 1)*(batch_q%length)/2) +                   &
            (n_mo - batch_q%length - batch_q%first + 1)*(batch_q%length)), &
            solver%n_cholesky)
!
      enddo
!
      call disk%close_file(solver%cholesky_ao_vectors)
      call disk%close_file(solver%cholesky_mo_vectors)
      call disk%close_file(cholesky_mo_vectors_seq, 'delete')
      call disk%close_file(solver%cholesky_ao_vectors_info, 'delete')
!
   end subroutine  construct_mo_cholesky_vecs
!
!
end module eri_cd_solver_class
!
!
