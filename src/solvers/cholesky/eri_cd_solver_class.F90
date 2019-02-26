module eri_cd_solver_class
!
!!
!!    Cholesky decomposition (CD) of electronic repulsion integrals (ERI) solver class
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
   use batching_index_class
   use molecular_system_class

!
   implicit none
!
!  Definition of the ERI-CD class 
!
   type :: eri_cd_solver
!
      character(len=100) :: tag = 'Cholesky decomposition of electronic repulsion integrals solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'

      character(len=500) :: description1 = 'Performs a Cholesky decomposition of the two-electron &
                                                  &electronic repulsion integrals in the atomic orbital basis,'

      character(len=500) :: description2 = '(ab|cd) = sum_J L_ab^J L_cd^J.'

      character(len=500) :: description3 = 'Once the Cholesky basis has been determined, the vectors &
                                                  &L^J are constructed and stored to disk. These may either be &
                                                  &used directly, or be transformed to the MO basis for use in &
                                                  &post-HF calculations.'
!
      real(dp) :: threshold         = 1.0D-8
      real(dp) :: span              = 1.0D-2
!
      integer :: max_qual      = 1000
!
      integer :: iteration     = 0
!
      logical :: one_center         = .false.
      logical :: construct_vectors  = .true.
!
      type(file) :: diagonal_info_target, diagonal_info_one_center, diagonal_info_construct
      type(file) :: cholesky_aux, cholesky_aux_inverse, cholesky_ao_vectors, cholesky_ao_vectors_info, cholesky_mo_vectors
      type(file) :: basis_shell_data
!
      integer :: n_cholesky
      integer :: n_sp_in_basis
      integer :: n_s, n_sp, n_ao, n_aop
!
      integer :: n_batches
!
   contains
!
      procedure :: prepare                                => prepare_eri_cd_solver
      procedure :: run                                    => run_eri_cd_solver
      procedure :: cleanup                                => cleanup_eri_cd_solver
!
      procedure :: invert_overlap_cholesky_vecs           => invert_overlap_cholesky_vecs_eri_cd_solver
      procedure :: cholesky_vecs_diagonal_test            => cholesky_vecs_diagonal_test_eri_cd_solver
      procedure :: diagonal_test                          => diagonal_test_eri_cd_solver
      procedure :: full_test_cholesky_vecs                => full_test_cholesky_vecs_cd_eri_solver
      procedure :: construct_significant_diagonal         => construct_significant_diagonal_eri_cd_solver
      procedure :: construct_significant_diagonal_atomic  => construct_significant_diagonal_atomic_eri_cd_solver
      procedure :: determine_auxilliary_cholesky_basis    => determine_auxilliary_cholesky_basis_eri_cd_solver
      procedure :: construct_overlap_cholesky_vecs        => construct_overlap_cholesky_vecs_eri_cd_solver
      procedure :: construct_cholesky_vectors             => construct_cholesky_vectors_eri_cd_solver
!
      procedure :: construct_mo_cholesky_vecs             => construct_mo_cholesky_vecs_cd_eri_solver
!
      procedure :: read_info                              => read_info_eri_cd_solver
      procedure :: print_banner                           => print_banner_eri_cd_solver
      procedure :: print_settings                         => print_settings_eri_cd_solver
!
      procedure :: construct_diagonal_batches             => construct_diagonal_batches_eri_cd_solver
      procedure :: construct_diagonal_from_batch_bases    => construct_diagonal_from_batch_bases_eri_cd_solver
      procedure :: append_bases                           => append_bases_eri_cd_solver
!
   end type eri_cd_solver
!
!
contains
!
!
   subroutine prepare_eri_cd_solver(solver, system)
!!
!!    Prepare ERI
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      type(molecular_system) :: system
!
      solver%n_batches = 1
!
      if (requested_section('cholesky')) then
!
         call solver%read_info()
!
      endif
!
      solver%n_aop   = system%get_n_aos()*(system%get_n_aos()+1)/2 ! Number of ao pairs packed
      solver%n_ao    = system%get_n_aos()
      solver%n_s     = system%get_n_shells()
      solver%n_sp    = solver%n_s*(solver%n_s + 1)/2              ! Number of shell pairs packed
!
!
      call solver%diagonal_info_target%init('target_diagonal', 'sequential', 'unformatted')
      call solver%diagonal_info_construct%init('construct_diagonal', 'sequential', 'unformatted')
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
   end subroutine prepare_eri_cd_solver
!
!
   subroutine run_eri_cd_solver(solver, system, screening_vector)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system) :: system
!
      real(dp), dimension(solver%n_ao), optional :: screening_vector
!
      real(dp):: s_determine_basis, e_determine_basis, s_build_vectors = 0, e_build_vectors = 0
      real(dp):: s_invert_time, e_invert_time, omp_get_wtime
!
      integer :: batch
!
      integer, dimension(:), allocatable :: n_cholesky_batches, n_sp_in_basis_batches
!
      type(file) :: batch_file_diag, batch_file_basis
!
      character(len=100) :: temp_name
!
      call solver%print_banner()
!
      call solver%print_settings()
!
      write(output%unit, '(/t6, a29, i13)') 'Total number of AOs:         ', system%get_n_aos()
      write(output%unit, '(t6, a29, i13)')  'Total number of shell pairs: ', solver%n_sp
      write(output%unit, '(t6, a29, i13)')  'Total number of AO pairs:    ', solver%n_aop
!
      write(output%unit, '(/t3, a38)') '- Preparing diagonal for decomposition'
!
      s_determine_basis = omp_get_wtime()
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
      else 
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
      endif
!
!
      if (solver%n_batches == 1) then
!
         call solver%determine_auxilliary_cholesky_basis(system, solver%diagonal_info_target, solver%basis_shell_data)
!
      else
!
         call mem%alloc(n_cholesky_batches, solver%n_batches)
         call mem%alloc(n_sp_in_basis_batches, solver%n_batches)
!
         n_cholesky_batches = 0
         n_sp_in_basis_batches = 0
!
         call solver%construct_diagonal_batches(system)
!
         write(output%unit, '(/t3, a31)') '- Decomposing batched diagonal:'
!
         do batch = 1, solver%n_batches 
!
            write(output%unit, '(/t3, a6, i3, a1)') 'Batch ', batch, ':'
!
            write(temp_name, '(a14, i4.4)')'diagonal_info_', batch
            call batch_file_diag%init(trim(temp_name), 'sequential', 'unformatted')

            write(temp_name, '(a11, i4.4)')'basis_info_', batch
            call batch_file_basis%init(trim(temp_name), 'sequential', 'unformatted')
!
            call solver%determine_auxilliary_cholesky_basis(system, batch_file_diag, batch_file_basis)
!
            n_cholesky_batches(batch)     = solver%n_cholesky
            n_sp_in_basis_batches(batch)  = solver%n_sp_in_basis
!
         enddo
!
         write(output%unit, '(/t3, a27)') '- Final decomposition step:'
!
         call solver%construct_diagonal_from_batch_bases(system, n_cholesky_batches, n_sp_in_basis_batches)
         call solver%determine_auxilliary_cholesky_basis(system, solver%diagonal_info_target, solver%basis_shell_data)
         !call solver%append_bases(system, n_cholesky_batches, n_sp_in_basis_batches)
!
         call mem%dealloc(n_cholesky_batches, solver%n_batches)
         call mem%dealloc(n_sp_in_basis_batches, solver%n_batches)
!
      endif
!
      e_determine_basis = omp_get_wtime()
      
      flush(output%unit)
!
      s_invert_time = omp_get_wtime()
!
      call solver%construct_overlap_cholesky_vecs(system)
      call solver%invert_overlap_cholesky_vecs()
!
      e_invert_time = omp_get_wtime()
               
      flush(output%unit)
!
      if (solver%construct_vectors) then
!
         s_build_vectors = omp_get_wtime()
!
         call solver%construct_cholesky_vectors(system)
!
         e_build_vectors = omp_get_wtime()
!
      endif
!
      write(output%unit, '(/t3, a)') '- Timings for Cholesky decomposition of electronic repulsion integrals:'
!
      write(output%unit, '(/t6, a53, f11.2)')'Wall time to determine auxiliary basis:              ', &
                                  e_determine_basis - s_determine_basis
      write(output%unit, '(t6, a53, f11.2)') 'Wall time to construct (J|K), decompose, and invert: ', &
                                  e_invert_time - s_invert_time
!
      if (solver%construct_vectors) &
         write(output%unit, '(t6, a53, f11.2)')'Wall time to build L_ab^J and test:                  ', &
                                  e_build_vectors - s_build_vectors                 
!
      flush(output%unit)
!
   end subroutine run_eri_cd_solver
!
!
   subroutine cleanup_eri_cd_solver(solver)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(solver%tag)
!
   end subroutine cleanup_eri_cd_solver
!
!
   subroutine construct_significant_diagonal_eri_cd_solver(solver, system, screening_vector)
!!
!!    Construct significant diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs two screened diagonals
!!
!!       1. Screened diagonal for decomposition D_αβ <= T   or   D_αβ * V_α * V_β <= T 
!!          The latter if optional screening vector is present
!!
!!       2. Screened diagonal for construction of cholesky vectros  sqrt( D_αβ * D_max ) <= T  
!!
!!    Writes all information to files target_diagonal and construct_diagonal
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      class(molecular_system) :: system
!
      real(dp), dimension(solver%n_ao), optional :: screening_vector
!
!     Local variables
!
      integer :: sp, n_sig_aop, n_sig_sp, current_sig_sp, n_construct_aop, n_construct_sp
!
      integer, dimension(:,:), allocatable :: sp_index, sig_sp_index
!
      integer, dimension(:), allocatable  :: ao_offsets
!
      real(dp), dimension(:), allocatable :: screening_vector_local, screening_vector_reduced, max_in_sp_diagonal
      real(dp), dimension(:,:), allocatable :: D_AB, D_AB_screen, construct_test
      real(dp), dimension(:), allocatable :: D_xy 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ABAB
!
      integer :: x, y, xy, xy_packed, A, B, I
!
      type(interval) :: A_interval, B_interval
!
      logical, dimension(:), allocatable :: sig_sp, construct_sp
!
      real(dp) :: max_diagonal
!
!     Prepare local screening vector
!
      call mem%alloc(screening_vector_local, solver%n_aop)
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
      call mem%alloc(sp_index, solver%n_sp, 2)
!
      sp = 0 ! Shell pair number
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
      call mem%alloc(max_in_sp_diagonal, solver%n_sp)
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, g_ABAB, D_AB, D_AB_screen) &
!$omp shared(sig_sp,  max_in_sp_diagonal) &
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
         call mem%alloc(g_ABAB, &
                  (A_interval%size),(B_interval%size), &
                  (A_interval%size),(B_interval%size))
!
         call system%ao_integrals%construct_ao_g_wxyz(g_ABAB, A, B, A, B)
!
         call mem%alloc(D_AB, (A_interval%size), (B_interval%size))
         call mem%alloc(D_AB_screen, (A_interval%size), (B_interval%size))
!
         do x = 1, (A_interval%size)
            do y = 1, (B_interval%size)
!
               D_AB_screen(x, y) = g_ABAB(x, y, x, y)&
                           *screening_vector_local(x + A_interval%first - 1)&
                           *screening_vector_local(y + B_interval%first - 1)
!
               D_AB(x, y) = g_ABAB(x, y, x, y)
!
            enddo
         enddo
!
         call mem%dealloc(g_ABAB, &
                  (A_interval%size), (B_interval%size), &
                  (A_interval%size), (B_interval%size))
!
!        Determine whether shell pair is significant
!
         sig_sp(I) = (is_significant(D_AB, (A_interval%size)*(B_interval%size), solver%threshold) .and. &
                      is_significant(D_AB_screen, (A_interval%size)*(B_interval%size), solver%threshold))
!
         max_in_sp_diagonal(I) = get_abs_max(D_AB, (A_interval%size)*(B_interval%size))
!
         call mem%dealloc(D_AB, (A_interval%size), (B_interval%size))
         call mem%dealloc(D_AB_screen, (A_interval%size), (B_interval%size))
!
      enddo
!$omp end parallel do
!
!     Prescreening for construction diagonal (Cauchy-Schwarz for final Cholesky vectors)
!
      max_diagonal = get_abs_max(max_in_sp_diagonal, solver%n_sp)
!
      call mem%dealloc(max_in_sp_diagonal, solver%n_sp)
!
      allocate(construct_sp(solver%n_sp))
      construct_sp = .false.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, g_ABAB, construct_test) &
!$omp shared(construct_sp) &
!$omp schedule(guided)
      do I = 1, solver%n_sp
!
         A = sp_index(I, 1)
         B = sp_index(I, 2)
!
         A_interval = system%shell_limits(A)
         B_interval = system%shell_limits(B)
!
!        Construct diagonal construct_test for the given shell pair
!
         call mem%alloc(g_ABAB, &
                  (A_interval%size), (B_interval%size), &
                  (A_interval%size), (B_interval%size))
!
         call system%ao_integrals%construct_ao_g_wxyz(g_ABAB, A, B, A, B)
!
         call mem%alloc(construct_test, (A_interval%size), (B_interval%size))
!
         do x = 1, (A_interval%size)
            do y = 1, (B_interval%size)
!
               construct_test(x, y) = sqrt(g_ABAB(x, y, x, y)*max_diagonal)
!
            enddo
         enddo
!
         call mem%dealloc(g_ABAB, &
                  (A_interval%size), (B_interval%size), &
                  (A_interval%size), (B_interval%size))
!
!           Determine whether shell pair should be constructed
!
            construct_sp(I) = is_significant(construct_test, (A_interval%size)*(B_interval%size), solver%threshold)
!
         call mem%dealloc(construct_test, (A_interval%size), (B_interval%size))
!
      enddo
!$omp end parallel do
!
!     Count number of screened AO pairs and shell pairs for both cases.
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
      n_construct_aop = 0 ! Number of AO pairs to construct
      n_construct_sp  = 0 ! Number of shell pairs to construct
!
      do I = 1, solver%n_sp
!
         if (construct_sp(I)) then
!
            A = sp_index(I, 1)
            B = sp_index(I, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            n_construct_aop = n_construct_aop + &
                           get_size_sp(A_interval, B_interval)
!
            n_construct_sp = n_construct_sp + 1
!
         endif
!
      enddo
!
      write(output%unit, '(/t6, a33, 2x, i11)')  'Significant shell pairs:         ', n_sig_sp
      write(output%unit, '(t6, a33, 2x, i11)')  'Significant AO pairs:            ', n_sig_aop
!
      write(output%unit, '(/t6, a33, 2x, i11)') 'Construct shell pairs:           ', n_construct_sp
      write(output%unit, '(t6, a33, 2x, i11)')  'Construct AO pairs:              ', n_construct_aop
      flush(output%unit)
!
!     Prepare for construction of diagonal and screening vector
!     
      call mem%alloc(ao_offsets, n_sig_sp)
!
      current_sig_sp = 0
!
      call mem%alloc(sig_sp_index, n_sig_sp, 2)
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
               ao_offsets(current_sig_sp + 1) = ao_offsets(current_sig_sp) + &
                           get_size_sp(A_interval, B_interval)
!
            endif
!
         endif
!
      enddo
!
      call mem%dealloc(sp_index, solver%n_sp, 2)
!
!     Construct significant diagonal and screening vector
!
      call mem%alloc(D_xy, n_sig_aop)
!
      call mem%alloc(screening_vector_reduced, n_sig_aop)
!
!     Note: allocated with length n_significant_sp + 1, last element is used for n_significant_aop
!     This is convenient because significant_sp_to_first_significant_aop will be used to calculate lengths.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, xy_packed, g_ABAB) &
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
         call mem%alloc(g_ABAB, &
               (A_interval%size), (B_interval%size), &
               (A_interval%size), (B_interval%size))
!
         call system%ao_integrals%construct_ao_g_wxyz(g_ABAB, A, B, A, B)
!
         if (A .eq. B) then
!
            do x = 1, A_interval%size
               do y = x, B_interval%size
!
                  xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                  D_xy(xy_packed + ao_offsets(I)) = g_ABAB(x, y, x, y)
                  screening_vector_reduced(xy_packed + ao_offsets(I)) = &
                                          screening_vector_local(x + A_interval%first - 1)*&
                                          screening_vector_local(y + B_interval%first - 1)
               enddo
            enddo
!
         else ! A ≠ B
!
            do x = 1, (A_interval%size)
               do y = 1, (B_interval%size)
!
                  xy = A_interval%size*(y - 1) + x
                  D_xy(xy + ao_offsets(I)) = g_ABAB(x, y, x, y)
                  screening_vector_reduced(xy + ao_offsets(I)) = &
                                      screening_vector_local(x + A_interval%first - 1)*&
                                      screening_vector_local(y + B_interval%first - 1)
!
               enddo
            enddo
!
         endif
!
         call mem%dealloc(g_ABAB, &
               (A_interval%size), (B_interval%size), &
               (A_interval%size), (B_interval%size))
!
      enddo
!$omp end parallel do   
!
      call mem%dealloc(screening_vector_local, solver%n_aop)
      call mem%dealloc(sig_sp_index, n_sig_sp, 2) 
      call mem%dealloc(ao_offsets, n_sig_sp)
!
!     Write info file for target diagonal containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!        4. Screening vector
!
      call disk%open_file(solver%diagonal_info_target, 'write', 'rewind')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
      write(solver%diagonal_info_target%unit) sig_sp
      write(solver%diagonal_info_target%unit) D_xy
      write(solver%diagonal_info_target%unit) screening_vector_reduced
!
!     Write info file for construct diagonal containing
!
!        1. number of shell pairs to construct, number of ao pairs to construct
!        2. construct_sp - vector of logicals to describe which shell pairs are to be constructed
!
      call disk%open_file(solver%diagonal_info_construct, 'write', 'rewind')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_construct%unit) n_construct_sp, n_construct_aop
      write(solver%diagonal_info_construct%unit) construct_sp
!
      call disk%close_file(solver%diagonal_info_target)
      call disk%close_file(solver%diagonal_info_construct)
!
      deallocate(sig_sp)
      call mem%dealloc(D_xy, n_sig_aop)
      call mem%dealloc(screening_vector_reduced, n_sig_aop)
!
   end subroutine construct_significant_diagonal_eri_cd_solver
!
!
   subroutine construct_significant_diagonal_atomic_eri_cd_solver(solver, system, screening_vector)
!!
!!    Construct significant diagonal atomic
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the significant diagonal for the given decomposition threshold within 
!!    the one-center approximation. Screening diagonal is optional argument for 
!!    additional screening.
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      class(molecular_system) :: system
!
      real(dp), dimension(solver%n_ao), optional :: screening_vector
!
      integer ::sp, n_sig_aop, n_sig_sp, current_sig_sp, n_construct_sp, n_construct_aop
!
      integer, dimension(:,:), allocatable :: sp_index, sig_sp_index
!
      integer, dimension(:), allocatable  :: ao_offsets
!
      real(dp), dimension(:), allocatable :: screening_vector_local, screening_vector_reduced, max_in_sp_diagonal
      real(dp), dimension(:), allocatable :: D_xy 
!
      real(dp), dimension(:,:), allocatable :: D_AB, D_AB_screen, construct_test
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ABAB
!
      integer :: x, y, xy, xy_packed, A, B, I
!
      type(interval) :: A_interval, B_interval
!
      logical, dimension(:), allocatable :: sig_sp, construct_sp
!
      real(dp) :: max_diagonal
!
      call mem%alloc(screening_vector_local, solver%n_aop)
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
      call mem%alloc(sp_index, solver%n_sp, 2)
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
      call mem%alloc(max_in_sp_diagonal, solver%n_sp)
!
!     Pre-screening of full diagonal
!
      allocate(sig_sp(solver%n_sp))
      sig_sp = .false.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, g_ABAB, D_AB, D_AB_screen) &
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
            call mem%alloc(g_ABAB, &
                     (A_interval%size), (B_interval%size), &
                     (A_interval%size), (B_interval%size))
!
            call system%ao_integrals%construct_ao_g_wxyz(g_ABAB, A, B, A, B)
!
            call mem%alloc(D_AB, (A_interval%size), (B_interval%size))
            call mem%alloc(D_AB_screen, (A_interval%size), (B_interval%size))
!
            do x = 1, (A_interval%size)
               do y = 1, (B_interval%size)
!
                  D_AB_screen(x, y) = g_ABAB(x, y, x, y)&
                              *screening_vector_local(x + A_interval%first - 1)&
                              *screening_vector_local(y + B_interval%first - 1)
!  
                  D_AB(x, y) = g_ABAB(x, y, x, y)
!
               enddo
            enddo
!
            call mem%dealloc(g_ABAB, &
                  (A_interval%size), (B_interval%size), &
                  (A_interval%size), (B_interval%size))
!
!           Determine whether shell pair is significant
!
            sig_sp(I) = (is_significant(D_AB, (A_interval%size)*(B_interval%size), solver%threshold) .and. &
                         is_significant(D_AB_screen, (A_interval%size)*(B_interval%size), solver%threshold))
!
            max_in_sp_diagonal(I) = get_abs_max(D_AB, (A_interval%size)*(B_interval%size))
!
            call mem%dealloc(D_AB, (A_interval%size), (B_interval%size))
            call mem%dealloc(D_AB_screen, (A_interval%size), (B_interval%size))
!
         endif
!
      enddo
!$omp end parallel do
!
!     Prescreening for construction diagonal (Cauchy-Schwarz for final Cholesky vectors)
!
      max_diagonal = get_abs_max(max_in_sp_diagonal, solver%n_sp)
!
      call mem%dealloc(max_in_sp_diagonal, solver%n_sp)
!
      allocate(construct_sp(solver%n_sp))
      construct_sp = .false.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, g_ABAB, construct_test) &
!$omp shared(construct_sp) &
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
         call mem%alloc(g_ABAB, &
                  (A_interval%size), (B_interval%size), &
                  (A_interval%size), (B_interval%size))
!
         call system%ao_integrals%construct_ao_g_wxyz(g_ABAB, A, B, A, B)
!
         call mem%alloc(construct_test, (A_interval%size), (B_interval%size))
!
         do x = 1, (A_interval%size)
            do y = 1, (B_interval%size)
!
               construct_test(x, y) = sqrt(g_ABAB(x, y, x, y)*max_diagonal)
!
            enddo
         enddo
!
         call mem%dealloc(g_ABAB, &
                  (A_interval%size), (B_interval%size), &
                  (A_interval%size), (B_interval%size))
!
!           Determine whether shell pair should be constructed
!
            construct_sp(I) = is_significant(construct_test, (A_interval%size)*(B_interval%size), solver%threshold)
!
         call mem%dealloc(construct_test, (A_interval%size), (B_interval%size))
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
      n_construct_aop = 0 ! Number of AO pairs to construct
      n_construct_sp  = 0 ! Number of shell pairs to construct
!
      do I = 1, solver%n_sp
!
         if (construct_sp(I)) then
!
            A = sp_index(I, 1)
            B = sp_index(I, 2)
!
            A_interval = system%shell_limits(A)
            B_interval = system%shell_limits(B)
!
            n_construct_aop = n_construct_aop + &
                           get_size_sp(A_interval, B_interval)
!
            n_construct_sp = n_construct_sp + 1
!
         endif
!
      enddo
!
      write(output%unit, '(/t6, a33, 2x, i11)')  'Significant shell pairs:         ', n_sig_sp
      write(output%unit, '(t6, a33, 2x, i11)')  'Significant AO pairs:            ', n_sig_aop
!
      write(output%unit, '(/t6, a33, 2x, i11)') 'Construct shell pairs:           ', n_construct_sp
      write(output%unit, '(t6, a33, 2x, i11)')  'Construct AO pairs:              ', n_construct_aop
      flush(output%unit)
!
!     Prepare for construction of diagonal and screening vector
!     
      call mem%alloc(ao_offsets, n_sig_sp)
      ao_offsets = 0
!
      current_sig_sp = 0
!
      call mem%alloc(sig_sp_index, n_sig_sp, 2)
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
               ao_offsets(current_sig_sp + 1) = ao_offsets(current_sig_sp) + &
                           get_size_sp(A_interval, B_interval)
!
            endif
!
         endif
!
      enddo
!
      call mem%dealloc(sp_index, solver%n_sp, 2)
!
!     Construct significant diagonal and screening vector
!
      call mem%alloc(D_xy, n_sig_aop)
!
      call mem%alloc(screening_vector_reduced, n_sig_aop)
!
!     Note: allocated with length n_significant_sp + 1, last element is used for n_significant_aop
!     This is convenient because significant_sp_to_first_significant_aop will be used to calculate lengths.
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, xy_packed, g_ABAB) &
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
         call mem%alloc(g_ABAB, &
               (A_interval%size), (B_interval%size), &
               (A_interval%size), (B_interval%size))
!
         call system%ao_integrals%construct_ao_g_wxyz(g_ABAB, A, B, A, B)
!
         if (A .eq. B) then
!
            do x = 1, A_interval%size
               do y = x, B_interval%size
!
                  xy = A_interval%size*(y - 1) + x
                  xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                  D_xy(xy_packed + ao_offsets(I)) = g_ABAB(x, y, x, y)
                  screening_vector_reduced(xy_packed + ao_offsets(I)) = &
                                          screening_vector_local(x + A_interval%first - 1)*&
                                          screening_vector_local(y + B_interval%first - 1)
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
                  D_xy(xy + ao_offsets(I)) = g_ABAB(x, y, x, y)
                  screening_vector_reduced(xy + ao_offsets(I)) = &
                                      screening_vector_local(x + A_interval%first - 1)*&
                                      screening_vector_local(y + B_interval%first - 1)
!
               enddo
            enddo
!
         endif
!
         call mem%dealloc(g_ABAB, &
               (A_interval%size), (B_interval%size), &
               (A_interval%size), (B_interval%size))
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(screening_vector_local, solver%n_aop)
      call mem%dealloc(sig_sp_index, n_sig_sp, 2)
      call mem%dealloc(ao_offsets, n_sig_sp)
!
!     Write info file for target diagonal containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!        4. Screening vector
!
      call disk%open_file(solver%diagonal_info_target, 'write', 'rewind')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
      write(solver%diagonal_info_target%unit) sig_sp
      write(solver%diagonal_info_target%unit) D_xy
      write(solver%diagonal_info_target%unit) screening_vector_reduced
!
!     Write info file for target diagonal containing
!
!        1. number of shell pairs to construct, number of ao pairs to construct
!        2. construct_sp - vector of logicals to describe which shell pairs are to be constructed
!
      call disk%open_file(solver%diagonal_info_construct, 'write', 'rewind')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_construct%unit) n_construct_sp, n_construct_aop
      write(solver%diagonal_info_construct%unit) construct_sp
!
      call disk%close_file(solver%diagonal_info_target)
      call disk%close_file(solver%diagonal_info_construct)
!
      deallocate(sig_sp)
      call mem%dealloc(D_xy, n_sig_aop)
      call mem%dealloc(screening_vector_reduced, n_sig_aop)
!
   end subroutine construct_significant_diagonal_atomic_eri_cd_solver
!
!
   subroutine construct_diagonal_batches_eri_cd_solver(solver, system)
!!
!!    Construct diagonal batches
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Divides the significant diagonal into batches and prepares for 
!!    decomposition
!!
      implicit none
!  
      class(eri_cd_solver) :: solver
!
      type(molecular_system) :: system
!
      integer :: n_sig_aop, n_sig_sp, n_sig_sp_batch, sp
!
      real(dp), dimension(:), allocatable :: D_xy, D_batch
!
      real(dp), dimension(:), allocatable :: screening_vector_batch, screening_vector
!
      logical, dimension(:), allocatable :: sig_sp, sig_sp_batch
!
      type(interval) :: A_interval, B_interval
!
      type(file) :: batch_file
!
      integer :: A, B, batch, batch_first, batch_last, batch_size, current_batch_size
      integer :: xy_first, xy_last
!
      character(len=100) :: temp_name
!
!     Read diagonal info file containing (name given as argument)
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call disk%open_file(solver%diagonal_info_target, 'read')
!
      rewind(solver%diagonal_info_target%unit)
!
      read(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
!
      call mem%alloc(D_xy, n_sig_aop)
      call mem%alloc(screening_vector, n_sig_aop)
      allocate(sig_sp(solver%n_sp))
!
      read(solver%diagonal_info_target%unit) sig_sp
      read(solver%diagonal_info_target%unit) D_xy
      read(solver%diagonal_info_target%unit) screening_vector
!
      call disk%close_file(solver%diagonal_info_target)
!
!     Calculate size of batches and the remainder
!
      batch_size = n_sig_aop/solver%n_batches
!
      allocate(sig_sp_batch(solver%n_sp))
!
      batch_first = 1
      batch_last = batch_size
!
      do batch = 1, solver%n_batches
!
!        Determine sig_sp_batch
!
         sig_sp_batch = .false.
!
         sp = 0        ! Shell pair number
         xy_first = 1
         xy_last = 0
         n_sig_sp_batch = 0
!
         do B = 1, solver%n_s
            do A = B, solver%n_s
!
               sp = sp + 1
!
               if (sig_sp(sp)) then 
!
                  A_interval = system%shell_limits(A)
                  B_interval = system%shell_limits(B)
!
                  xy_last = xy_last + get_size_sp(A_interval, B_interval)
!
                  if ((xy_last .ge. batch_first) .and. (xy_first .le. batch_last)) then
!
                     sig_sp_batch(sp) = .true.
                     n_sig_sp_batch = n_sig_sp_batch + 1
!
                     if (xy_last .gt. batch_last) then 
!
                        batch_last = xy_last
!
                     endif
!
                  endif
!
                  xy_first = xy_first + get_size_sp(A_interval, B_interval)
!
               endif
!
            enddo
         enddo     
!
         current_batch_size = batch_last - batch_first + 1
!
         call mem%alloc(D_batch, current_batch_size)
         call mem%alloc(screening_vector_batch, current_batch_size)
!
         D_batch(:) = D_xy(batch_first : batch_last)
!
         screening_vector_batch(:) = screening_vector_batch(batch_first : batch_last)
!
!        Write info file for batch diagonal containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!        4. Screening vector
!
         write(temp_name, '(a14, i4.4)')'diagonal_info_', batch
         call batch_file%init(trim(temp_name), 'sequential', 'unformatted')
!
         call disk%open_file(batch_file, 'write', 'rewind')
         rewind(batch_file%unit)
!
         write(output%unit, '(/t6, a40, i3, a1)')         'Significant AO and shell pairs in batch ', batch, ':'
!
         write(output%unit, '(/t9, a33, 2x, i11)')  'Significant shell pairs:         ', n_sig_sp_batch
         write(output%unit, '(t9, a33, 2x, i11)')  'Significant AO pairs:            ', current_batch_size
!
         flush(output%unit)
!
!
         write(batch_file%unit) n_sig_sp_batch, current_batch_size
         write(batch_file%unit) sig_sp_batch
         write(batch_file%unit) D_batch
         write(batch_file%unit) screening_vector_batch
!
         call disk%close_file(batch_file)
!
         call mem%dealloc(D_batch, current_batch_size)
         call mem%dealloc(screening_vector_batch, current_batch_size)
!
         batch_first = batch_last + 1  
         batch_last  = batch_size*(batch + 1)
!
         if ((batch + 1) == solver%n_batches) batch_last = n_sig_aop
!
         if (batch_last .lt. batch_first) call output%error_msg('batch size is too small.')
!
      enddo
!
      deallocate(sig_sp_batch)
      deallocate(sig_sp)
      call mem%dealloc(D_xy, n_sig_aop)
      call mem%dealloc(screening_vector, n_sig_aop)
!
   end subroutine construct_diagonal_batches_eri_cd_solver
!
!
   subroutine construct_diagonal_from_batch_bases_eri_cd_solver(solver, system, n_cholesky_batches, n_sp_in_basis_batches)
!!
!!    Construct diagonal from batch bases
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
!!    Constructs the final diagonal from the bases obtained from diagonal batches.
!!    Called as preparation for final decomposition step.
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system) :: system
!
      integer, dimension(solver%n_batches), intent(in) :: n_cholesky_batches
      integer, dimension(solver%n_batches), intent(in) :: n_sp_in_basis_batches
!
      integer :: n_cholesky_total, n_sp_in_basis_total, J, I, n_sig_aop, n_sig_sp
      integer :: A_shell, B_shell, sp, alpha_in_A, beta_in_B, alpha_beta_in_AB, aop, batch
      integer :: n_basis_aop_in_AB_total, n_basis_aop_in_AB_offset, current_offset, current_offset_old
      integer :: count_sig, n_cholesky_offset, n_sig_aop_old, n_sig_sp_old, n_sp_in_basis_offset
!
      type(file) :: batch_file
!
      integer, dimension(:), allocatable :: alpha, beta, alpha_beta, sorted_alpha, sorted_beta, sorted_alpha_beta
      integer, dimension(:), allocatable :: index_alpha_beta, alpha_beta_offset, alpha_beta_offset_old
      integer, dimension(:), allocatable :: A, B, AB, sorted_A, sorted_B, sorted_AB
      integer, dimension(:), allocatable :: sorted_n_basis_aop_in_AB, n_basis_aop_in_AB, index_AB
!
      integer, dimension(:,:), allocatable :: basis_shell_info, cholesky_basis
!
      logical, dimension(:), allocatable :: sig_sp, sig_sp_old
!
      type(interval) :: A_interval, B_interval
!
      real(dp), dimension(:), allocatable :: D, D_old, screening_vector, screening_vector_old
!
      character(len=100) :: temp_name
!
      n_cholesky_total = 0
      n_sp_in_basis_total = 0
!
      do batch = 1, solver%n_batches
!
         n_cholesky_total    = n_cholesky_total + n_cholesky_batches(batch) 
         n_sp_in_basis_total = n_sp_in_basis_total + n_sp_in_basis_batches(batch) 
!
      enddo
!
!     Read and paste together basis information 
!     from the different batches 
!
      call mem%alloc(alpha, n_cholesky_total)
      call mem%alloc(beta, n_cholesky_total)
      call mem%alloc(alpha_beta, n_cholesky_total)
!
      call mem%alloc(A, n_sp_in_basis_total)
      call mem%alloc(B, n_sp_in_basis_total)
      call mem%alloc(AB, n_sp_in_basis_total)
      call mem%alloc(n_basis_aop_in_AB, n_sp_in_basis_total)
!
      n_sp_in_basis_offset = 0
      n_cholesky_offset = 0
!
      do batch = 1, solver%n_batches
!
!        Basis_shell_data file order:
!  
!           1. number of sps in basis 
!           2. basis_shell_info
!           3. cholesky_basis
         
         write(temp_name, '(a11, i4.4)') 'basis_info_', batch
         call batch_file%init(trim(temp_name), 'sequential', 'unformatted')
!  
         call disk%open_file(batch_file, 'read')
         rewind(batch_file%unit)
!
         call mem%alloc(basis_shell_info, n_sp_in_basis_batches(batch), 4)
         call mem%alloc(cholesky_basis, n_cholesky_batches(batch), 3)
!  
         read(batch_file%unit) 
         read(batch_file%unit) basis_shell_info
         read(batch_file%unit) cholesky_basis
!  
         call disk%close_file(batch_file)
!
         do J = 1, n_cholesky_batches(batch)
!
            alpha(n_cholesky_offset + J)      = cholesky_basis(J, 1)
            beta(n_cholesky_offset + J)       = cholesky_basis(J, 2)
            alpha_beta(n_cholesky_offset + J) = cholesky_basis(J, 3)
!
         enddo
!
         do sp = 1, n_sp_in_basis_batches(batch)
!
            A(n_sp_in_basis_offset + sp)                 = basis_shell_info(sp, 1)
            B(n_sp_in_basis_offset + sp)                 = basis_shell_info(sp, 2)
            AB(n_sp_in_basis_offset + sp)                = basis_shell_info(sp, 3)
            n_basis_aop_in_AB(n_sp_in_basis_offset + sp) = basis_shell_info(sp, 4)
!
         enddo
!
         n_sp_in_basis_offset = n_sp_in_basis_offset + n_sp_in_basis_batches(batch)
         n_cholesky_offset    = n_cholesky_offset + n_cholesky_batches(batch)
!
         call mem%dealloc(basis_shell_info, n_sp_in_basis_batches(batch), 4)
         call mem%dealloc(cholesky_basis, n_cholesky_batches(batch), 3)
!
      enddo
!
!     Sort the arrays according to an alphabeta and an AB ordering 
!     from smallest to largest
!
      call mem%alloc(index_AB, n_sp_in_basis_total)
      call quicksort_with_index_ascending_int(AB, index_AB, n_sp_in_basis_total)
!
      call mem%alloc(index_alpha_beta, n_cholesky_total)
      call quicksort_with_index_ascending_int(alpha_beta, index_alpha_beta, n_cholesky_total)
!
      call mem%alloc(sorted_alpha, n_cholesky_total)
      call mem%alloc(sorted_beta, n_cholesky_total)
      call mem%alloc(sorted_alpha_beta, n_cholesky_total)
!
      sorted_alpha_beta = alpha_beta
      call mem%dealloc(alpha_beta, n_cholesky_total)
!
      do J = 1, n_cholesky_total
!
         sorted_alpha(J) = alpha(index_alpha_beta(J))
         sorted_beta(J)  = beta(index_alpha_beta(J))
!
      enddo
!
      call mem%dealloc(alpha, n_cholesky_total)
      call mem%dealloc(beta, n_cholesky_total)
!
      call mem%alloc(sorted_A, n_sp_in_basis_total)
      call mem%alloc(sorted_B, n_sp_in_basis_total)
      call mem%alloc(sorted_AB, n_sp_in_basis_total)
      call mem%alloc(sorted_n_basis_aop_in_AB, n_sp_in_basis_total)
!
      sorted_AB = AB
!
      call mem%dealloc(AB, n_sp_in_basis_total)
!
      do J = 1, n_sp_in_basis_total
!
         sorted_A(J) = A(index_AB(J))
         sorted_B(J) = B(index_AB(J))
         sorted_n_basis_aop_in_AB(J) = n_basis_aop_in_AB(index_AB(J))
!
      enddo
!
      call mem%dealloc(A, n_sp_in_basis_total)
      call mem%dealloc(B, n_sp_in_basis_total)
      call mem%dealloc(n_basis_aop_in_AB, n_sp_in_basis_total)
!
      call mem%dealloc(index_alpha_beta, n_cholesky_total)
      call mem%dealloc(index_AB, n_sp_in_basis_total)
!
!     Construct significant shell pair logical array,
!     and count the number of significant AO and shell pairs 
!
      allocate(sig_sp(solver%n_sp))
!
      sig_sp = .false.
!
      n_sig_sp = 0
      n_sig_aop = 0
!
      I = 1
      sp = 0
!    
      do B_shell = 1, solver%n_s
         do A_shell = B_shell, solver%n_s
!
            sp = sp + 1
!
            if (sp == sorted_AB(I)) then 
!
               A_interval = system%shell_limits(A_shell)
               B_interval = system%shell_limits(B_shell)
!
               sig_sp(sp) = .true.
               n_sig_sp = n_sig_sp + 1
               n_sig_aop = n_sig_aop + get_size_sp(A_interval, B_interval)
!
               I = I + 1
!
            endif
!
         enddo
      enddo
!
!     Read old diagonal from file, along with old sig_sp logical array, related info.,
!     and screening vector (old refers here to the initially screened diagonal)
!
      call mem%alloc(D, n_sig_aop)
      call mem%alloc(screening_vector, n_sig_aop)
!
      call disk%open_file(solver%diagonal_info_target, 'read')
!
      rewind(solver%diagonal_info_target%unit)
!
      read(solver%diagonal_info_target%unit) n_sig_sp_old, n_sig_aop_old
!
      call mem%alloc(D_old, n_sig_aop_old)
      call mem%alloc(screening_vector_old, n_sig_aop_old)
      allocate(sig_sp_old(solver%n_sp))
!
      read(solver%diagonal_info_target%unit) sig_sp_old
      read(solver%diagonal_info_target%unit) D_old
      read(solver%diagonal_info_target%unit) screening_vector_old
!
      call disk%close_file(solver%diagonal_info_target)
!
!     Copy the correct elements of the initial D into the new D using cholesky basis array.
!     We precalculate alpha beta offsets both old and new, then copy afterwards.
!
      call mem%alloc(alpha_beta_offset, n_sig_sp)
      call mem%alloc(alpha_beta_offset_old, n_sig_sp)
!
      count_sig = 0
      current_offset = 0
      current_offset_old = 0
!
      sp = 0
!
      do B_shell = 1, solver%n_s
         do A_shell = B_shell, solver%n_s
!
            sp = sp + 1
!
            if (sig_sp_old(sp)) then
!
               A_interval = system%shell_limits(A_shell)
               B_interval = system%shell_limits(B_shell)
!
               if (sig_sp(sp)) then 
!
                  count_sig = count_sig + 1
!
                  alpha_beta_offset_old(count_sig) = current_offset_old
                  alpha_beta_offset(count_sig)     = current_offset
!
                  current_offset = current_offset + get_size_sp(A_interval, B_interval)
!
               endif
!
               current_offset_old = current_offset_old + get_size_sp(A_interval, B_interval)
!
            endif
!
         enddo
!
      enddo
!
      I = 0
      count_sig = 0
      n_basis_aop_in_AB_offset = 0
!
      do while (I .lt. n_sp_in_basis_total)
!
         I = I + 1
         count_sig = count_sig + 1
!
         n_basis_aop_in_AB_total = sorted_n_basis_aop_in_AB(I)
!
         A_interval = system%shell_limits(sorted_A(I))
         B_interval = system%shell_limits(sorted_B(I))
!
         do aop = 1, n_basis_aop_in_AB_total 
!
            alpha_in_A = sorted_alpha(aop + n_basis_aop_in_AB_offset) - A_interval%first + 1
            beta_in_B  = sorted_beta(aop + n_basis_aop_in_AB_offset) - B_interval%first + 1
!
            if (sorted_A(I) == sorted_B(I)) then
!
               alpha_beta_in_AB = max(alpha_in_A, beta_in_B)*(max(alpha_in_A, beta_in_B) - 3)/2 + alpha_in_A + beta_in_B
!
            else
!
               alpha_beta_in_AB = A_interval%size*(beta_in_B - 1) + alpha_in_A
!
            endif
!
            D(alpha_beta_in_AB + alpha_beta_offset(count_sig)) = &
                                    D_old(alpha_beta_in_AB + alpha_beta_offset_old(count_sig))
!        
            screening_vector(alpha_beta_in_AB + alpha_beta_offset(count_sig)) = &
                                    screening_vector_old(alpha_beta_in_AB + alpha_beta_offset_old(count_sig))
!
         enddo
!
         n_basis_aop_in_AB_offset = n_basis_aop_in_AB_offset + n_basis_aop_in_AB_total 
!
      enddo
!
      call mem%dealloc(alpha_beta_offset, n_sig_sp)
      call mem%dealloc(alpha_beta_offset_old, n_sig_sp)
!
      call mem%dealloc(D_old, n_sig_aop_old)
      call mem%dealloc(screening_vector_old, n_sig_aop_old)
!
      call mem%dealloc(sorted_alpha, n_cholesky_total)
      call mem%dealloc(sorted_beta, n_cholesky_total)
      call mem%dealloc(sorted_alpha_beta, n_cholesky_total)
!
      call mem%dealloc(sorted_A, n_sp_in_basis_total)
      call mem%dealloc(sorted_B, n_sp_in_basis_total)
      call mem%dealloc(sorted_AB, n_sp_in_basis_total)
      call mem%dealloc(sorted_n_basis_aop_in_AB, n_sp_in_basis_total)
!
      deallocate(sig_sp_old)
!
!     Write info file for target diagonal containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!        4. Screening vector
!
      write(output%unit, '(/t6, a)')  'Significant AO and shell pairs in final decomposition:'
!
      write(output%unit, '(/t6, a33, 2x, i11)')  'Significant shell pairs:         ', n_sig_sp
      write(output%unit, '(t6, a33, 2x, i11)')  'Significant AO pairs:            ', n_sig_aop
!
      flush(output%unit)
!
      call disk%open_file(solver%diagonal_info_target, 'write', 'rewind')
      rewind(solver%diagonal_info_target%unit)
!
      write(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
      write(solver%diagonal_info_target%unit) sig_sp
      write(solver%diagonal_info_target%unit) D
      write(solver%diagonal_info_target%unit) screening_vector
!
      call disk%close_file(solver%diagonal_info_target)
!
      deallocate(sig_sp)
!
      call mem%dealloc(D, n_sig_aop)
      call mem%dealloc(screening_vector, n_sig_aop)
!
   end subroutine construct_diagonal_from_batch_bases_eri_cd_solver
!
!
   subroutine append_bases_eri_cd_solver(solver, n_cholesky_batches, n_sp_in_basis_batches)
!!
!!    Append bases
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Appends bases from different diagonal batches, to be used
!!    if system routine is used to decompose directly the bases from the batches
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      integer, dimension(solver%n_batches), intent(in) :: n_cholesky_batches
      integer, dimension(solver%n_batches), intent(in) :: n_sp_in_basis_batches
!
      integer :: n_cholesky_total, n_sp_in_basis_total, J
      integer :: sp, batch
      integer :: n_cholesky_offset, n_sp_in_basis_offset
!
      type(file) :: batch_file
!
      integer, dimension(:,:), allocatable :: basis_shell_info, cholesky_basis
      integer, dimension(:,:), allocatable :: cholesky_full, basis_shell_info_full
!
      character(len=100) :: temp_name
!
      n_cholesky_total = 0
      n_sp_in_basis_total = 0
!
      do batch = 1, solver%n_batches
!
         n_cholesky_total = n_cholesky_total + n_cholesky_batches(batch) 
         n_sp_in_basis_total = n_sp_in_basis_total + n_sp_in_basis_batches(batch) 
!
      enddo
!
      call mem%alloc(cholesky_full, n_cholesky_total, 3)
      call mem%alloc(basis_shell_info_full, n_sp_in_basis_total, 4)
!
      n_sp_in_basis_offset = 0
      n_cholesky_offset = 0
!
      do batch = 1, solver%n_batches
!
!        Read basis_shell_data file containing
!  
!           1. number shell pairs in basis
!           2. basis_shell_info
!           3. cholesky_basis
         
         write(temp_name, '(a11, i4.4)')'basis_info_', batch
         call batch_file%init(trim(temp_name), 'sequential', 'unformatted')
!  
         call disk%open_file(batch_file, 'read')
         rewind(batch_file%unit)
!
         call mem%alloc(basis_shell_info, n_sp_in_basis_batches(batch), 4)
         call mem%alloc(cholesky_basis, n_cholesky_batches(batch), 3)
!  
         read(batch_file%unit) 
         read(batch_file%unit) basis_shell_info
         read(batch_file%unit) cholesky_basis
!  
         call disk%close_file(batch_file)
!
         do J = 1, n_cholesky_batches(batch)
!
            cholesky_full(n_cholesky_offset + J, :) = cholesky_basis(J, :)
!
         enddo
!
         do sp = 1, n_sp_in_basis_batches(batch)
!
            basis_shell_info_full(n_sp_in_basis_offset + sp, :) = basis_shell_info(sp, :)
!
         enddo
!
         n_sp_in_basis_offset = n_sp_in_basis_offset + n_sp_in_basis_batches(batch)
         n_cholesky_offset = n_cholesky_offset + n_cholesky_batches(batch)
!
         call mem%dealloc(basis_shell_info, n_sp_in_basis_batches(batch), 4)
         call mem%dealloc(cholesky_basis, n_cholesky_batches(batch), 3)
!
      enddo
!
      call disk%open_file(solver%basis_shell_data, 'write', 'rewind')
!
      write(solver%basis_shell_data%unit) n_sp_in_basis_total
!
      write(solver%basis_shell_data%unit) basis_shell_info_full
      write(solver%basis_shell_data%unit) cholesky_full
!
      solver%n_cholesky = n_cholesky_total
!
      call disk%close_file(solver%basis_shell_data)
!
      call mem%dealloc(cholesky_full, n_cholesky_total, 3)
      call mem%dealloc(basis_shell_info_full, n_sp_in_basis_total, 4)
!
   end subroutine append_bases_eri_cd_solver
!
!
   subroutine determine_auxilliary_cholesky_basis_eri_cd_solver(solver, system, diagonal_info, basis_info)
!!
!!    Determine auxiliary cholesky basis
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the elements of the auxiliary basis for an RI-type expansion
!!    by Cholesky decomposition
!!
!!
      implicit none
!
      class(eri_cd_solver), intent(inout) :: solver
!
      class(molecular_system), intent(in) :: system
!
      type(file), intent(in) :: diagonal_info, basis_info
!
!     Local variables
!
!     Integers
!
      integer :: n_sig_sp, n_sig_aop
      integer :: n_new_sig_sp, n_new_sig_aop, current_new_sig_sp
      integer :: n_qual_sp, n_qual_aop, n_previous_qual_aop, n_qual_aop_in_sp
      integer :: sp, current_sig_sp
      integer :: first_sig_aop, last_sig_aop, aop
      integer :: A, B, AB, AB_sp
      integer :: C, D, CD_sp
      integer :: I, J
      integer :: w, x, y, z
      integer :: xy, xy_packed, xy_max, wx, wx_packed, yz
      integer :: sig_neg
      integer :: first, last
      integer :: first_x, first_y
      integer :: current_qual, qual
      integer :: n_cholesky_in_node, n_new_cholesky
      integer :: n_sp_in_basis
      integer :: sig_sp_counter
      integer :: sp_in_basis
!
!     Integer allocatable arrays
!
      integer, dimension(:), allocatable :: sig_sp_to_first_sig_aop         ! Maps significant shell pair to first ao pair
      integer, dimension(:), allocatable :: new_sig_sp_to_first_sig_aop     ! Maps significant shell pair to first ao pair
      integer, dimension(:), allocatable :: sorted_max_sig_sp               ! Index array for sorting shell pairs according to their maximum values
      integer, dimension(:), allocatable :: sorted_qual_aop_in_sp_indices   ! Index array for sorting the qualified ao pairs in shell pair
      integer, dimension(:), allocatable :: n_qual_aop_in_prev_sps          ! Offsets for omp-loop, number of qualified ao pairs in preceding shell pair
      integer, dimension(:), allocatable :: qual_max                        ! Index list containing order in which qualified diagonals are selected in decomposition
      integer, dimension(:), allocatable :: sig_sp_to_previous_sig_sp       ! Maps significant shell pair indices to significant shell pair indices of last iteration, used for reduction
!      
      integer, dimension(:,:), allocatable :: sig_sp_to_shells              ! Maps significant shell pair to shells
      integer, dimension(:,:), allocatable :: new_sig_sp_to_shells          ! Maps significant shell pair to shells   
      integer, dimension(:,:), allocatable :: sig_aop_to_aos                ! Maps significant ao pair to aos
      integer, dimension(:,:), allocatable :: new_sig_aop_to_aos            ! Maps significant ao pair to aos      
      integer, dimension(:,:), allocatable :: qual_sp                       ! List of qualified shell pairs
      integer, dimension(:,:), allocatable :: qual_sp_copy                  ! List of qualified shell pairs, copy used to reduce size
      integer, dimension(:,:), allocatable :: qual_aop                      ! List of qualified ao pairs
      integer, dimension(:,:), allocatable :: qual_aop_copy                 ! List of qualified ao pairs, copy used to reduce size
      integer, dimension(:,:), allocatable :: cholesky_basis                ! ao and ao pair indices of the elements of the cholesky basis
      integer, dimension(:,:), allocatable :: cholesky_basis_new            ! ao and ao pair indices of the elements of the cholesky basis, written to file at end of routine
      integer, dimension(:,:), allocatable :: basis_shell_info_full         ! Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: basis_shell_info              ! Info on shells containing elements of the basis, written to file at end of routine
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
      real(dp), dimension(:), allocatable :: D_xy                             ! Array for eri diagonal elements
      real(dp), dimension(:), allocatable :: D_xy_new                         ! Array for eri diagonal elements, used for reduction
      real(dp), dimension(:), allocatable :: approx_diagonal_accumulative     ! Array for accumulating approximate diagonal
      real(dp), dimension(:), allocatable :: max_in_sig_sp                    ! Maximum in each significant shell pair
      real(dp), dimension(:), allocatable :: sorted_qual_aop_in_sp            ! Sorted qualified ao pair in shell pair
      real(dp), dimension(:), allocatable :: screening_vector                 ! Screening vector for diagonal
      real(dp), dimension(:), allocatable :: screening_vector_new             ! Screening vector for diagonal, used for reduction
!      
      real(dp), dimension(:,:), allocatable :: g_wxyz                         ! Array for eri
      real(dp), dimension(:,:), allocatable :: cholesky_tmp                   ! Array used for dgemm, reordered copy of cholesky vectors of current batch of qualified
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ABCD                     ! Array for eri for shell pairs AB and CD
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
      call mem%alloc(D_xy, n_sig_aop)
      call mem%alloc(screening_vector, n_sig_aop)
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
      call mem%alloc(sig_sp_to_first_sig_aop, n_sig_sp + 1)
      sig_sp_to_first_sig_aop = 0
!
      sig_sp_to_first_sig_aop(n_sig_sp + 1) = n_sig_aop + 1
!
      call mem%alloc(sig_sp_to_shells, n_sig_sp, 2) ! [A, B]
      sig_sp_to_shells = 0
!
      call mem%alloc(sig_aop_to_aos, n_sig_aop, 2) ! [alpha, beta]
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
               sig_sp_to_first_sig_aop(current_sig_sp) = first_sig_aop
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
      write(output%unit, '(/t3, a)')&
      'Iter.  #Sign. ao pairs / shell pairs   Max diagonal    #Qualified    #Cholesky    Cholesky array size'
      write(output%unit, '(t3, a)') &
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
         call mem%alloc(max_in_sig_sp, n_sig_sp)
!
         max_in_sig_sp = zero
!
         do sp = 1, n_sig_sp
!
!           Get first and last indices of shell pair
!
            first = sig_sp_to_first_sig_aop(sp)
            last  = sig_sp_to_first_sig_aop(sp + 1) - 1
!
!           Determine the largest elements
!
            do I = first, last
!
               if (D_xy(I) .gt. max_in_sig_sp(sp)) then
!
                  max_in_sig_sp(sp) = D_xy(I)
!
               endif
!
            enddo
!
         enddo
!
!        Sort from largest to smallest and determine an index array of sorting
!
         call mem%alloc(sorted_max_sig_sp, n_sig_sp)
         sorted_max_sig_sp = 0
!
         call quicksort_with_index_descending(max_in_sig_sp, sorted_max_sig_sp, n_sig_sp)
!
         D_max_full  = max_in_sig_sp(1)
         n_qual_aop  = 0
         n_qual_sp   = 0
!
         call mem%dealloc(max_in_sig_sp, n_sig_sp)
!
         call mem%alloc(qual_aop, solver%max_qual, 3)
         call mem%alloc(qual_sp, solver%n_sp, 3)
         qual_sp = 0
         qual_aop = 0
!
         do sp = 1, n_sig_sp
!
            current_sig_sp = sorted_max_sig_sp(sp)
!
            first_sig_aop = sig_sp_to_first_sig_aop(current_sig_sp)
            last_sig_aop  = sig_sp_to_first_sig_aop(current_sig_sp + 1) - 1
!
            n_qual_aop_in_sp = 0
!
            do aop = first_sig_aop, last_sig_aop
!
               if ((D_xy(aop) .ge. solver%span*D_max_full) .and. (n_qual_aop .lt. solver%max_qual)) then
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
               call mem%alloc(sorted_qual_aop_in_sp_indices, n_qual_aop_in_sp)
               call mem%alloc(sorted_qual_aop_in_sp, n_qual_aop_in_sp)
!
               call get_n_highest(n_qual_aop_in_sp, last_sig_aop - first_sig_aop + 1, &
                                 D_xy(first_sig_aop:last_sig_aop), sorted_qual_aop_in_sp, &
                                 sorted_qual_aop_in_sp_indices)
!
               n_previous_qual_aop = (n_qual_aop - n_qual_aop_in_sp)
!
               do aop = 1, n_qual_aop_in_sp
!
                  qual_aop(aop + n_previous_qual_aop, 1) = sig_aop_to_aos(sorted_qual_aop_in_sp_indices(aop) &
                                                               + first_sig_aop - 1, 1)
!
                  qual_aop(aop + n_previous_qual_aop, 2) = sig_aop_to_aos(sorted_qual_aop_in_sp_indices(aop) &
                                                               + first_sig_aop - 1, 2)
!
                  qual_aop(aop + n_previous_qual_aop, 3) = sorted_qual_aop_in_sp_indices(aop) &
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
               call mem%dealloc(sorted_qual_aop_in_sp_indices, n_qual_aop_in_sp)
               call mem%dealloc(sorted_qual_aop_in_sp, n_qual_aop_in_sp)
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
         call mem%dealloc(sorted_max_sig_sp, n_sig_sp)
!
!        Cut out the qualified parts of the aop and sp lists
!
         call mem%alloc(qual_aop_copy, n_qual_aop, 3)
         call mem%alloc(qual_sp_copy, n_qual_sp, 3)
!
         qual_aop_copy(:, :) = qual_aop(1 : n_qual_aop, :)
         qual_sp_copy(:, :)  = qual_sp(1 : n_qual_sp, :)
!
         call mem%dealloc(qual_aop, solver%max_qual, 3)
         call mem%dealloc(qual_sp, solver%n_sp, 3)
!
         call mem%alloc(qual_aop, n_qual_aop, 3)
         call mem%alloc(qual_sp, n_qual_sp, 3)
!
         qual_aop    = qual_aop_copy
         qual_sp     = qual_sp_copy
!
         call mem%dealloc(qual_aop_copy, n_qual_aop, 3)
         call mem%dealloc(qual_sp_copy, n_qual_sp, 3)
!
!        Prepare to construct g_wxyz in parallelized loop
!
         call mem%alloc(n_qual_aop_in_prev_sps, n_qual_sp)
         n_qual_aop_in_prev_sps = 0
!
         do CD_sp = 1, n_qual_sp - 1
!
             n_qual_aop_in_prev_sps(CD_sp + 1) = n_qual_aop_in_prev_sps(CD_sp) + qual_sp(CD_sp, 3)
!
         enddo
!
!        Construct g_wxyz
!
         call mem%alloc(g_wxyz, n_sig_aop, n_qual_aop)
!
!$omp parallel do &
!$omp private(AB_sp, CD_sp, A, B, A_interval, B_interval, C, D, C_interval, D_interval, &
!$omp  aop, w, x, y, z, wx, yz, wx_packed, g_ABCD, n_qual_aop_in_sp) &
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
               call mem%alloc(g_ABCD, &
                              (A_interval%size), (B_interval%size), &
                              (C_interval%size), (D_interval%size))
!
               call system%ao_integrals%construct_ao_g_wxyz(g_ABCD, A, B, C, D)
!
                 do aop = 1, n_qual_aop_in_sp
!
                    y = qual_aop(aop + n_qual_aop_in_prev_sps(CD_sp), 1)
                    z = qual_aop(aop + n_qual_aop_in_prev_sps(CD_sp), 2)
!
                    if (A == B) then
!
                       do w = 1, A_interval%size
                          do x = w, B_interval%size
!
                             wx_packed = (max(w,x)*(max(w,x)-3)/2) + w + x
!
                             g_wxyz(sig_sp_to_first_sig_aop(AB_sp) + wx_packed - 1, aop + n_qual_aop_in_prev_sps(CD_sp)) &
                                   = g_ABCD(w, x, y - C_interval%first + 1, z - D_interval%first + 1)
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
                             g_wxyz(sig_sp_to_first_sig_aop(AB_sp) + wx - 1, aop + n_qual_aop_in_prev_sps(CD_sp)) &
                                   = g_ABCD(w, x, y - C_interval%first + 1, z - D_interval%first + 1)
!
                          enddo
                       enddo
!
                    endif
!
                 enddo
!
                  call mem%dealloc(g_ABCD, &
                                    (A_interval%size), (B_interval%size), &
                                    (C_interval%size), (D_interval%size))
!
            enddo
!
         enddo ! cd_sp
!$omp end parallel do
!
         call mem%dealloc(n_qual_aop_in_prev_sps, n_qual_sp)
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
            call mem%alloc(cholesky_basis, solver%n_cholesky + n_qual_aop, 3)
            cholesky_basis(1 : solver%n_cholesky, :) = cholesky_basis_new(:, :)
            call mem%dealloc(cholesky_basis_new, solver%n_cholesky, 3)
!
         else
!
            call mem%alloc(cholesky_basis, n_qual_aop, 3)
            cholesky_basis = 0
!
         endif
!
         call mem%alloc(approx_diagonal_accumulative, n_sig_aop)
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
         call mem%alloc(qual_max, n_qual_aop)
!
         do while ((current_qual .lt. n_qual_aop) .and. construct_more_choleskys)
!
            current_qual = current_qual + 1
!
            D_max = zero
            xy_max = 0
!
            do qual = 1, n_qual_aop
!
               xy = qual_aop(qual, 3)
!
               if (D_xy(xy) - approx_diagonal_accumulative(xy) .gt. D_max) then
!
                  qual_max(current_qual) = qual
                  D_max    = D_xy(xy) - approx_diagonal_accumulative(xy)
                  xy_max = xy
!
               endif
!
            enddo
!
            if ((D_max*screening_vector(xy_max) .gt. solver%threshold) .and. &
               (D_max .gt. solver%threshold)  .and. &
               (D_max .ge. solver%span*D_max_full)) then
!
               cholesky_basis(solver%n_cholesky + current_qual, 1) = qual_aop(qual_max(current_qual), 1)
               cholesky_basis(solver%n_cholesky + current_qual, 2) = qual_aop(qual_max(current_qual), 2)
!
               A = system%basis2shell(qual_aop(qual_max(current_qual), 1))
               B = system%basis2shell(qual_aop(qual_max(current_qual), 2))
!
               cholesky_basis(solver%n_cholesky + current_qual, 3) = get_sp_from_shells(A, B, solver%n_s)
!
               cholesky_new(: , current_qual) = g_wxyz(:, qual_max(current_qual))
!
               if (current_qual .gt. 1) then
!
                  call mem%alloc(cholesky_tmp, 1, current_qual - 1)
!
                  cholesky_tmp(1, :) = cholesky_new(qual_aop(qual_max(current_qual), 3), 1 : current_qual - 1)
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
                  if (D_xy(xy) == zero) then
!
                     cholesky_new(xy, current_qual) = zero
!
                  endif
!
               enddo
!
               do xy = 1, n_sig_aop
!
                  approx_diagonal_accumulative(xy) = approx_diagonal_accumulative(xy) + cholesky_new(xy, current_qual)**2
!
               enddo
!
               D_xy(qual_aop(qual_max(current_qual), 3)) = zero
               approx_diagonal_accumulative(qual_aop(qual_max(current_qual), 3))  = zero
!
               do xy = 1, n_sig_aop
!
                  if (D_xy(xy) - approx_diagonal_accumulative(xy) .lt. zero) then
!
                     if (abs(D_xy(xy) - approx_diagonal_accumulative(xy)) .gt. 1.0d-10) then
                        if (write_warning) then
                           write(output%unit, '(t6, a)') 'Warning: Found significant negative diagonal! '
                           write_warning = .false.
                        endif
                        sig_neg = sig_neg + 1
                     endif
!
                     D_xy(xy) = zero
                     approx_diagonal_accumulative(xy) = zero
!
                  elseif ((D_xy(xy) - approx_diagonal_accumulative(xy))*screening_vector(xy) .lt. solver%threshold .or. &
                     (D_xy(xy) - approx_diagonal_accumulative(xy)).lt. solver%threshold) then
!
                     D_xy(xy) = zero
                     approx_diagonal_accumulative(xy) = zero
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
            D_xy(xy) = D_xy(xy) - approx_diagonal_accumulative(xy)
!
         enddo
!
         call mem%dealloc(approx_diagonal_accumulative, n_sig_aop)
!
         call mem%dealloc(qual_max, n_qual_aop)
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
               first = sig_sp_to_first_sig_aop(sig_sp_counter)
               last  = sig_sp_to_first_sig_aop(sig_sp_counter + 1) - 1
!
               new_sig_sp(sig_sp_counter) = (is_significant(D_xy(first:last), &
                                                last - first + 1, solver%threshold, &
                                                screening_vector(first:last) ) .and. &
                                             is_significant(D_xy(first:last), &
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
            call mem%alloc(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1)
            new_sig_sp_to_first_sig_aop = 0
!
            call mem%alloc(sig_sp_to_previous_sig_sp, n_sig_sp + 1) ! 1 2 3 4 ... n_sig_sp, n_sig_sp + 1
            sig_sp_to_previous_sig_sp(n_sig_sp + 1) = n_sig_sp + 1
!
            current_new_sig_sp    = 1
            n_new_sig_aop = 0
            first_sig_aop = 1
!
            do sp = 1, n_sig_sp
!
               sig_sp_to_previous_sig_sp(sp) = sp
!
               if (new_sig_sp(sp)) then
!
                  A = sig_sp_to_shells(sp, 1)
                  B = sig_sp_to_shells(sp, 2)
!
                  A_interval = system%shell_limits(A)
                  B_interval = system%shell_limits(B)
!
                  new_sig_sp_to_first_sig_aop(current_new_sig_sp) = first_sig_aop
!
                  first_sig_aop = first_sig_aop + get_size_sp(A_interval, B_interval)
                  n_new_sig_aop = first_sig_aop - 1
!
                  current_new_sig_sp    = current_new_sig_sp + 1

               endif
!
            enddo
!
            new_sig_sp_to_first_sig_aop(current_new_sig_sp) = n_new_sig_aop + 1
!
            call mem%alloc(new_sig_aop_to_aos, n_new_sig_aop, 2)
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
            call mem%alloc(new_sig_sp_to_shells, n_new_sig_sp, 2)
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
            call mem%dealloc(sig_sp_to_previous_sig_sp, n_sig_sp + 1)
            call mem%dealloc(sig_sp_to_shells, n_sig_sp, 2)
            call mem%alloc(sig_sp_to_shells, n_new_sig_sp, 2)
!
            sig_sp_to_shells = new_sig_sp_to_shells
            call mem%dealloc(new_sig_sp_to_shells, n_new_sig_sp, 2)
!
            call mem%alloc(D_xy_new, n_new_sig_aop)
!
           call reduce_vector(D_xy,                     &
                             D_xy_new,                  &
                             sig_sp_to_first_sig_aop,   &
                             new_sig_sp,                &
                             n_sig_sp,                  &
                             n_sig_aop,                 &
                             n_new_sig_aop)
!
            call mem%dealloc(D_xy, n_sig_aop)
            call mem%alloc(D_xy, n_new_sig_aop)
!
            call dcopy(n_new_sig_aop, D_xy_new, 1, D_xy, 1)
!
            call mem%dealloc(D_xy_new, n_new_sig_aop)
!
            call mem%alloc(screening_vector_new, n_new_sig_aop)
!
           call reduce_vector(screening_vector,         &
                             screening_vector_new,      &
                             sig_sp_to_first_sig_aop,   &
                             new_sig_sp,                &
                             n_sig_sp,                  &
                             n_sig_aop,                 &
                             n_new_sig_aop)
!
            call mem%dealloc(screening_vector, n_sig_aop)
            call mem%alloc(screening_vector, n_new_sig_aop)
!
            call dcopy(n_new_sig_aop, screening_vector_new, 1, screening_vector, 1)
!
            call mem%dealloc(screening_vector_new, n_new_sig_aop)
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
            call mem%alloc(cholesky_basis_new, solver%n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : solver%n_cholesky + n_new_cholesky, :)
            call mem%dealloc(cholesky_basis, solver%n_cholesky + n_qual_aop, 3)
!
   !        Deallocate old lists & reallocate + copy over new lists
!
            deallocate(new_sig_sp)
!
            call mem%dealloc(sig_sp_to_first_sig_aop, n_sig_sp + 1)
            call mem%alloc(sig_sp_to_first_sig_aop, n_new_sig_sp + 1)
            sig_sp_to_first_sig_aop = new_sig_sp_to_first_sig_aop
            call mem%dealloc(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1)
!
            call mem%dealloc(sig_aop_to_aos, n_sig_aop, 2)
            call mem%alloc(sig_aop_to_aos, n_new_sig_aop, 2)
            sig_aop_to_aos = new_sig_aop_to_aos
            call mem%dealloc(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            n_sig_sp = n_new_sig_sp
            n_sig_aop = n_new_sig_aop
!
            solver%n_cholesky = solver%n_cholesky + n_new_cholesky
!
            call mem%dealloc(qual_aop, n_qual_aop, 3)
            call mem%dealloc(qual_sp, n_qual_sp, 3)
!
            write(output%unit, '(t3, i4, 8x, i9, 1x, a1, i9, 6x, e12.5, 4x, i4, 8x, i7, 8x, i13)') &
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
            call mem%dealloc(D_xy, n_sig_aop)
!
            call mem%alloc(cholesky_basis_new, solver%n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : solver%n_cholesky + n_new_cholesky, :)
            call mem%dealloc(cholesky_basis, solver%n_cholesky + n_qual_aop, 3)
            call mem%dealloc(screening_vector, n_sig_aop)
!
            solver%n_cholesky = solver%n_cholesky + n_new_cholesky
!
            done = .true.
!
            write(output%unit, '(t3, i4, 8x, i9, 1x, a1, i9, 6x, e12.5, 4x, i4, 8x, i7, 8x, i13)') &
            solver%iteration, 0,'/',0, D_max_full, n_qual_aop, solver%n_cholesky, 0
            flush(output%unit)
!
         endif
!
         call cpu_time(e_reduce_time)
         full_reduce_time = full_reduce_time + e_reduce_time - s_reduce_time
!
      enddo ! while not done
      write(output%unit, '(t3, a)') &
      '-------------------------------------------------------------------------------------------------------'
      flush(output%unit)
!
!     Timings
!
      call cpu_time(e_select_basis_time)
!
      if (sig_neg .gt. 0) write(output%unit,'(/t6, a42, i7)')'Number of significant negative diagonals: ', sig_neg
!
!     Prepare info on basis
!
!     Construct a list of all shell pairs (and shells) that contain elements of the basis
!     and how many elements of the basis they contain
!
      call mem%alloc(basis_shell_info_full, solver%n_sp, 4) ! A, B, AB, n_basis_aops_in_sp
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
      call mem%alloc(basis_shell_info, n_sp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_sp_in_basis, :)
      call mem%dealloc(basis_shell_info_full, solver%n_sp, 4)
!
!     Write basis_shell_data file containing
!
!        1. number shell pairs in basis
!        2. basis_shell_info
!        3. cholesky_basis
!
      call disk%open_file(basis_info, 'write', 'rewind')
!
      write(basis_info%unit) n_sp_in_basis
!
      write(basis_info%unit) basis_shell_info
      write(basis_info%unit) cholesky_basis_new
!
      solver%n_sp_in_basis = n_sp_in_basis
!
      call disk%close_file(basis_info)
!
      call mem%dealloc(basis_shell_info, n_sp_in_basis, 4)
      call mem%dealloc(cholesky_basis_new, solver%n_cholesky, 3)
!
   end subroutine determine_auxilliary_cholesky_basis_eri_cd_solver
!
!
   subroutine construct_overlap_cholesky_vecs_eri_cd_solver(solver, system)
!!
!!    Construct overlap cholesky vectors
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the overlap matrix (J|K) of the auxiliary basis 
!!    and cholesky decomposes it.
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
      integer :: n_sp_in_basis, sp_in_basis
      integer :: n_vectors
      integer :: current_aop_in_sp
      integer :: A, B, C, D, AB, AB_sp, CD_sp
      integer :: I, J, K, L, KL
      integer :: w, x, y, z, wx, yz
!
!     Integer allocatable arrays
!
      integer, dimension(:,:), allocatable :: basis_shell_info        ! Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: basis_shell_info_full   ! Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: cholesky_basis          ! ao and ao pair indices of the elements of the cholesky basis
      integer, dimension(:,:), allocatable :: cholesky_basis_updated  ! ao and ao pair indices of the elements of the cholesky basis
      integer, dimension(:,:), allocatable :: basis_aops_in_CD_sp     ! basis ao pairs in shell pair CD
      integer, dimension(:,:), allocatable :: basis_aops_in_AB_sp     ! basis ao pairs in shell pair AB
!
      integer, dimension(:), allocatable :: keep_vectors
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
      real(dp), dimension(:,:), allocatable :: integrals_auxiliary
!
      real(dp), dimension(:), allocatable :: work  ! work array for LAPACK
!
      integer :: info
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
      call mem%alloc(basis_shell_info, n_sp_in_basis, 4)
      call mem%alloc(cholesky_basis, solver%n_cholesky, 3)
!
      read(solver%basis_shell_data%unit) basis_shell_info
      read(solver%basis_shell_data%unit) cholesky_basis
!
      call disk%close_file(solver%basis_shell_data, 'delete')
!
!     Construct integrals (J | J')
!
      call mem%alloc(integrals_auxiliary, solver%n_cholesky, solver%n_cholesky)
!
!$omp parallel do &
!$omp private(AB_sp, CD_sp, A, B, A_interval, B_interval, C, D, C_interval, D_interval, &
!$omp w, x, y, z, wx, yz, g_AB_CD, I, J, K, L, KL,&
!$omp current_aop_in_sp, basis_aops_in_CD_sp, basis_aops_in_AB_sp) &
!$omp shared(integrals_auxiliary, cholesky_basis, basis_shell_info) &
!$omp schedule(guided)
      do CD_sp = 1, n_sp_in_basis
!
         C = basis_shell_info(CD_sp, 1)
         D = basis_shell_info(CD_sp, 2)
!
         C_interval = system%shell_limits(C)
         D_interval = system%shell_limits(D)
!
         call mem%alloc(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
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
            call mem%alloc(basis_aops_in_AB_sp, basis_shell_info(AB_sp, 4), 3)
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
            call system%ao_integrals%construct_ao_g_wxyz(g_AB_CD, A, B, C, D)
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
!
                  integrals_auxiliary(K, L) = g_AB_CD(wx, yz)
                  integrals_auxiliary(L, K) = g_AB_CD(wx, yz)
!
               enddo
            enddo
!
            call mem%dealloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
            call mem%dealloc(basis_aops_in_AB_sp, basis_shell_info(AB_sp, 4), 3)
!
         enddo ! AB
!
         call mem%dealloc(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
      enddo ! CD
!$omp end parallel do
!
      call mem%dealloc(basis_shell_info, n_sp_in_basis, 4)
!
      n_vectors = 0
      allocate(keep_vectors(solver%n_cholesky))
!
      call cpu_time(s_decomp_time)
!
      allocate(work(2*solver%n_cholesky))
!
!     DPSTRF computes the Cholesky factorization with complete pivoting
!     of a real symmetric positive semidefinite matrix.
!
      call dpstrf('L',                 &
            solver%n_cholesky,         &
            integrals_auxiliary,       &
            solver%n_cholesky,         &
            keep_vectors,              &
            n_vectors,                 &
            solver%threshold*1.0d-1,   &
            work,                      &
            info)
!
      deallocate(work)
!
      call cpu_time(e_decomp_time)
!
      write(output%unit, '(/t6, a)') 'Done decomposing (J|K)!'
      flush(output%unit)
!
      call mem%alloc(cholesky_basis_updated, n_vectors, 3)
!
!$omp parallel do private(I)
      do I = 1, n_vectors
!
         cholesky_basis_updated(I, :) = cholesky_basis(keep_vectors(I), :)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(cholesky_basis, solver%n_cholesky, 3)
      deallocate(keep_vectors)
!
!     Write cholesky_aux file containing
!
!        1. Cholesky vectors L_JK
!
      call disk%open_file(solver%cholesky_aux, 'write', 'rewind')
!
      write(solver%cholesky_aux%unit) integrals_auxiliary(1:n_vectors, 1:n_vectors)
!
      call disk%close_file(solver%cholesky_aux)
!
      call mem%dealloc(integrals_auxiliary, solver%n_cholesky, solver%n_cholesky)
!
      solver%n_cholesky = n_vectors
!
      write(output%unit, '(/t3, a, i10)')'- Summary of Cholesky decomposition of electronic repulsion integrals: '
      write(output%unit, '(/t6, a, i10)')'Final number of Cholesky vectors: ', solver%n_cholesky
      flush(output%unit)
!
!     Update the basis_shell_info array which contains information of which shell pairs (and shells)
!     contain elements of the basis and how many elements of the basis they contain.
!
      call mem%alloc(basis_shell_info_full, solver%n_sp, 4) ! A, B, AB, n_basis_aops_in_sp
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
      call mem%alloc(basis_shell_info, n_sp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_sp_in_basis, :)
      call mem%dealloc(basis_shell_info_full, solver%n_sp, 4)
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
      call mem%dealloc(basis_shell_info, n_sp_in_basis, 4)
      call mem%dealloc(cholesky_basis_updated, n_vectors, 3)
!
      call cpu_time(e_build_basis_time)
!
   end subroutine construct_overlap_cholesky_vecs_eri_cd_solver
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
      real(dp):: s_invert_time, e_invert_time
!
      real(dp), dimension(:,:), allocatable :: cholesky_inverse
!
      integer :: I, info
!
      call cpu_time(s_invert_time)
!
!     Read Cholesky vectors of auxiliary basis overlap
!
      call disk%open_file(solver%cholesky_aux, 'read')
      rewind(solver%cholesky_aux%unit)
!
      call mem%alloc(cholesky_inverse, solver%n_cholesky, solver%n_cholesky)
!
      read(solver%cholesky_aux%unit) cholesky_inverse
!
      call disk%close_file(solver%cholesky_aux, 'delete')
!
!     Invert cholesky vectors
!
      call DTRTRI('l','n', solver%n_cholesky, cholesky_inverse, solver%n_cholesky, info)
!
      if (info /= 0) call output%error_msg('Error: matrix inversion failed!', info)
!
      write(output%unit, '(/t6, a)') 'Done inverting L_JK!'
      flush(output%unit)
!
!     Write inverse Cholesky vectors of auxiliary basis overlap
!
!        1. n_cholesky: Number of elements of the basis
!        2. cholesky_inverse (n_cholesky, n_cholesky)
!
      call disk%open_file(solver%cholesky_aux_inverse, 'write', 'rewind')
!
!     Write out columns, but only the parts that are not rubbish left over from dpstrf
!
      do I = 1, solver%n_cholesky
!
         write(solver%cholesky_aux_inverse%unit) cholesky_inverse(1 + (I - 1) : solver%n_cholesky, I)
!
      enddo
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
   subroutine construct_cholesky_vectors_eri_cd_solver(solver, system)
!!
!!    Construct Cholesky vectors
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the Cholesky vectors L_xy_J = sum_K L_K_J^-1 (K | xy)
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
      integer :: A, B, AB_sp, C, D, CD_sp
      integer :: w, x, y, z, wx, yz, wx_packed
      integer :: L, J, I
      integer :: n_construct_sp, n_construct_aop
      integer :: n_sp_in_basis, last_sp_included, sp_counter
      integer :: current_aop_in_sp
      integer :: n_AB_included, n_AB_included_current
      integer :: rec_offset
      integer :: size_AB, size_AB_current
!
!     Integer allocatable arrays
!
      integer, dimension(:,:), allocatable :: basis_shell_info     ! Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: AB_info              ! Info on offsets and shells for OMP-loop [offset, A, B]
      integer, dimension(:,:), allocatable :: basis_aops_in_CD_sp  ! Basis ao pairs in shell pair CD
      integer, dimension(:,:), allocatable :: cholesky_basis       ! Info on cholesky basis
!
!     Reals
!
      real(dp) :: s_construct_time, e_construct_time, full_construct_time
      real(dp) :: s_build_vectors_time
!
!     Real allocatable arrays
!
      real(dp), dimension(:,:), allocatable :: g_wx_L, g_AB_CD
      real(dp), dimension(:,:), allocatable :: L_wx_J
      real(dp), dimension(:,:), allocatable :: aux_chol_inverse, aux_chol_inverse_transpose
!
!     Logicals
!
      logical :: done, found_size
!
!     Logical allocatable arrays
!
      logical, dimension(:), allocatable :: construct_sp
!
!     Intervals
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
      call cpu_time(s_build_vectors_time)
!
!     Read diagonal info
!
      call disk%open_file(solver%diagonal_info_construct, 'read')
      rewind(solver%diagonal_info_construct%unit)
!
      allocate(construct_sp(solver%n_sp))
!
      read(solver%diagonal_info_construct%unit) n_construct_sp, n_construct_aop
      read(solver%diagonal_info_construct%unit) construct_sp
!
      call disk%close_file(solver%diagonal_info_construct)
!
!     Read inverse of cholesky vectors of auxiliary overlap, L_JK^-1
!
      call disk%open_file(solver%cholesky_aux_inverse, 'read')
      rewind(solver%cholesky_aux_inverse%unit)
!
      call mem%alloc(aux_chol_inverse, solver%n_cholesky, solver%n_cholesky)
      aux_chol_inverse = zero
!
      do I = 1, solver%n_cholesky
!
         read(solver%cholesky_aux_inverse%unit) aux_chol_inverse(1 + (I - 1) : solver%n_cholesky, I)
!
      enddo
!
      call disk%close_file(solver%cholesky_aux_inverse)
!
!     Transpose L_JK^-1
!
      call mem%alloc(aux_chol_inverse_transpose, solver%n_cholesky, solver%n_cholesky)
!
      call trans(aux_chol_inverse, aux_chol_inverse_transpose, solver%n_cholesky)
!
      call mem%dealloc(aux_chol_inverse, solver%n_cholesky, solver%n_cholesky)
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
         last_sp_included        = 0
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
               if (construct_sp(sp_counter)) then
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
         if (.not. found_size) then ! No batching
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
         rewind(solver%basis_shell_data%unit)
!
         read(solver%basis_shell_data%unit) n_sp_in_basis
!
         call mem%alloc(basis_shell_info, n_sp_in_basis, 4)
         call mem%alloc(cholesky_basis, solver%n_cholesky, 3)
!
         read(solver%basis_shell_data%unit) basis_shell_info
         read(solver%basis_shell_data%unit) cholesky_basis
!
         call disk%close_file(solver%basis_shell_data)
!
         call mem%alloc(AB_info, n_AB_included, 3) ! [offset, A, B]
         AB_info = zero
!
         sp_counter = 0
         AB_sp = 0
!
         do B = 1, solver%n_s
            do A = B, solver%n_s
!
               AB_sp = AB_sp + 1
!
               if (construct_sp(AB_sp) .and. AB_sp .le. last_sp_included) then
!
                  sp_counter = sp_counter + 1
                  construct_sp(AB_sp) = .false.       ! So that this shell pair is not selected again
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
!        Construct g_J_yz = (J | yz)
!
         call mem%alloc(g_wx_L, size_AB, solver%n_cholesky)
!
!$omp parallel do &
!$omp private(AB_sp, CD_sp, I, A, B, A_interval, &
!$omp B_interval, C, D, C_interval, D_interval, &
!$omp basis_aops_in_CD_sp, current_aop_in_sp, g_AB_CD, &
!$omp w, x, y, z, wx, yz, wx_packed, L, J) &
!$omp shared(g_wx_L, AB_info, basis_shell_info, cholesky_basis) &
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
               call mem%alloc(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
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
              call mem%alloc(g_AB_CD, &
                       (A_interval%size)*(B_interval%size), &
                       (C_interval%size)*(D_interval%size))
!
              g_AB_CD = zero
!
              call system%ao_integrals%construct_ao_g_wxyz(g_AB_CD, A, B, C, D)
!
               if (A == B) then
!
                  do w = 1, A_interval%size
                     do x = w, B_interval%size
!
                        wx_packed = (max(w,x)*(max(w,x)-3)/2) + w + x 
                        wx = A_interval%size*(x-1) + w
!
                        do J = 1, basis_shell_info(CD_sp, 4)
                           y = basis_aops_in_CD_sp(J, 1)
                           z = basis_aops_in_CD_sp(J, 2)
                           L = basis_aops_in_CD_sp(J, 3)
                           yz = C_interval%size*(z-1)+y

                           g_wx_L(wx_packed + AB_info(AB_sp, 1), L) = g_AB_CD(wx, yz)
!
                           enddo
                        enddo
                     enddo
!
                  else
!
                     do w = 1, A_interval%size
                        do x = 1, B_interval%size
                           do J = 1, basis_shell_info(CD_sp, 4)
                              y = basis_aops_in_CD_sp(J, 1)
                              z = basis_aops_in_CD_sp(J, 2)
                              L = basis_aops_in_CD_sp(J, 3)
!
                              yz = C_interval%size*(z-1) + y
!
                              wx = A_interval%size*(x-1) + w
!
                              g_wx_L(wx + AB_info(AB_sp, 1), L) = g_AB_CD(wx, yz)
!
                           enddo
                        enddo
                     enddo
!
                  endif
!
                  call mem%dealloc(g_AB_CD,                 &
                     (A_interval%size)*(B_interval%size),   &
                     (C_interval%size)*(D_interval%size))
!
               call mem%dealloc(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
            enddo ! CD
!
         enddo ! AB
!$omp end parallel do
!
         call mem%dealloc(AB_info, n_AB_included, 3)
!
         call cpu_time(s_construct_time)
!
!        L_wx_J = sum_L (wx | L) (L | J)^-T 
!
         call mem%alloc(L_wx_J, size_AB, solver%n_cholesky)
!
         call dgemm('N', 'N',                      &
                       size_AB,                    &
                       solver%n_cholesky,          &
                       solver%n_cholesky,          &
                       one,                        &
                       g_wx_L,                     & 
                       size_AB,                    &
                       aux_chol_inverse_transpose, & 
                       solver%n_cholesky,          &
                       zero,                       &
                       L_wx_J,                     &
                       size_AB)
!
         call mem%dealloc(g_wx_L, size_AB, solver%n_cholesky)
!
         call cpu_time(e_construct_time)
         full_construct_time = full_construct_time + e_construct_time - s_construct_time
!
!        Write vectors to file
!
         do J = 1, solver%n_cholesky
!
            write(solver%cholesky_ao_vectors%unit) (L_wx_J(I, J), I = 1, size_AB)
!
         enddo
!
         rec_offset = rec_offset + size_AB
!
         call mem%dealloc(L_wx_J, size_AB, solver%n_cholesky)
!
         done = .true.
!
         do I = 1, solver%n_sp
            if (construct_sp(I)) then
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
      call mem%dealloc(aux_chol_inverse_transpose, solver%n_cholesky, solver%n_cholesky)
      call mem%dealloc(cholesky_basis, solver%n_cholesky, 3)
      call mem%dealloc(basis_shell_info, n_sp_in_basis, 4)
      deallocate(construct_sp)
!
   end subroutine construct_cholesky_vectors_eri_cd_solver
!
!
   subroutine diagonal_test_eri_cd_solver(solver, system)
!!
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
      type(molecular_system) :: system
!
      real(dp), dimension(:,:), allocatable :: D_red, D, g_AB_AB
!
      real(dp) :: max_diff
!
      integer :: n_sig_sp, n_sig_aop, AB_offset, x, y, xy, xy_red, xy_full, sp, A, B
!
      type(interval) :: A_interval, B_interval
!
      logical, dimension(:), allocatable :: sig_sp
!
!     Read diagonal information
!
      call disk%open_file(solver%diagonal_info_target, 'read')
      rewind(solver%diagonal_info_target%unit)
!
      read(solver%diagonal_info_target%unit) n_sig_sp, n_sig_aop
!
      call mem%alloc(D_red, n_sig_aop, 1)
!
      allocate(sig_sp(solver%n_sp))
!
      read(solver%diagonal_info_target%unit) sig_sp
      read(solver%diagonal_info_target%unit) D_red
!
      call disk%close_file(solver%diagonal_info_target)
!
      call mem%alloc(D, solver%n_ao**2, 1)
      D = zero
!
      AB_offset = 0
      sp = 0
!
      do B = 1, solver%n_s 
!
         B_interval = system%shell_limits(B)
!
         do A = B, solver%n_s
!
            A_interval = system%shell_limits(A)
!
            sp = sp + 1
!
            call mem%alloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
               g_AB_AB = zero
               call system%ao_integrals%construct_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
               do x = 1, (A_interval%size)
                  do y = 1, (B_interval%size)
!
                     xy = (A_interval%size)*(y-1)+x
                     xy_full = solver%n_ao*(y + B_interval%first - 2) + x + A_interval%first - 1
!
                     D(xy_full, 1) = g_AB_AB(xy, xy)
!
                  enddo
               enddo
!
               call mem%dealloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
            if (sig_sp(sp)) then
!
               do x = 1, (A_interval%size)
                  do y = 1, (B_interval%size)
!
                     xy_full = solver%n_ao*(y + B_interval%first - 2) + x + A_interval%first - 1
!
                     if (A .ne. B) then 
!
                        xy_red = (A_interval%size)*(y-1)+x
!
                     else
!
                        xy_red  = (max(y,x)*(max(y,x)-3)/2) + y + x
!
                     endif
!
                     D(xy_full, 1) = D(xy_full, 1) - D_red(xy_red + AB_offset, 1)
!
                  enddo
               enddo
!
               AB_offset = AB_offset + get_size_sp(A_interval, B_interval)
!
            endif
         enddo
      enddo
!
!     Calculate maximal difference and minimal difference
!
      max_diff = get_abs_max(D, solver%n_ao**2)
!
      call mem%dealloc(D, solver%n_ao**2, 1)
!
      write(output%unit, '(/t2, a)')'- Testing the Screened diagonal:'
!
      write(output%unit, '(/t6, a, e12.4)')'Maximal difference between screened and actual diagonal: ', max_diff
      flush(output%unit)
!
   end subroutine diagonal_test_eri_cd_solver
!
!
   subroutine cholesky_vecs_diagonal_test_eri_cd_solver(solver, system)
!!
!!    Cholesky vectors diagonal test
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Tests the decomposition by 
!!       1. finding the largest element of (D_sig - D_approx)
!!       2. finding the smallest (largest negative) element of (D_sig - D_approx)
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system), intent(in) :: system
!
      real(dp), dimension(:,:), allocatable :: D_xy, L_K_yz, g_AB_AB
!
      real(dp) :: ddot, max_diff, min_diff
!
      integer :: aop, n_construct_sp, n_construct_aop, J, I, size_AB, AB_offset, A, B, current_construct_sp
      integer :: x, y, xy, xy_packed, sp
!
      integer, dimension(:,:), allocatable :: construct_sp_index, ao_offsets
!
      character(len=40) :: line
!
      logical, dimension(:), allocatable :: construct_sp
!
      type(interval) :: A_interval, B_interval
!
!     Read diagonal information
!
      call disk%open_file(solver%diagonal_info_construct, 'read')
!
      rewind(solver%diagonal_info_construct%unit)
!
      read(solver%diagonal_info_construct%unit) n_construct_sp, n_construct_aop
!
      allocate(construct_sp(solver%n_sp))
!
      read(solver%diagonal_info_construct%unit) construct_sp
!
      call disk%close_file(solver%diagonal_info_construct)
!
!     Prepare for construction of diagonal
!
      sp = 0        ! Shell pair number
!    
      call mem%alloc(ao_offsets, n_construct_sp, 1)
      ao_offsets = 0
!
      current_construct_sp = 0
!
      call mem%alloc(construct_sp_index, n_construct_sp, 2)
!
      do B = 1, solver%n_s
         do A = B, solver%n_s
!
            sp = sp + 1
!
            if (construct_sp(sp)) then
!
               current_construct_sp = current_construct_sp + 1
!
               A_interval = system%shell_limits(A)
               B_interval = system%shell_limits(B)
!     
               construct_sp_index(current_construct_sp, 1) = A
               construct_sp_index(current_construct_sp, 2) = B
!
               if (current_construct_sp .lt. n_construct_sp) then
!
                  ao_offsets(current_construct_sp + 1, 1) = ao_offsets(current_construct_sp, 1) + &
                           get_size_sp(A_interval, B_interval)
!
               endif
!
            endif
!
         enddo
      enddo
!
!     Construct significant diagonal and screening vector
!
      call mem%alloc(D_xy, n_construct_aop, 1)
      D_xy = zero
!
!$omp parallel do &
!$omp private(I, A, B, A_interval, B_interval, x, y, xy, xy_packed, g_AB_AB) &
!$omp shared(D_xy, ao_offsets) &
!$omp schedule(guided)
      do I = 1, n_construct_sp
!
         A = construct_sp_index(I, 1)
         B = construct_sp_index(I, 2)
!
         A_interval = system%shell_limits(A)
         B_interval = system%shell_limits(B)
!
         call mem%alloc(g_AB_AB, &
               (A_interval%size)*(B_interval%size), &
               (A_interval%size)*(B_interval%size))
!
         g_AB_AB = zero
         call system%ao_integrals%construct_ao_g_wxyz(g_AB_AB, A, B, A, B)
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
      call mem%dealloc(construct_sp_index, n_construct_sp, 2)
      call mem%dealloc(ao_offsets, n_construct_sp, 1)
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
            D_xy(I + AB_offset, 1) = D_xy(I + AB_offset, 1) - ddot(solver%n_cholesky, L_K_yz(1, I), 1, L_K_yz(1, I), 1)
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
      max_diff = get_abs_max(D_xy, n_construct_aop)
!
      min_diff =1.0d5
!
      do aop = 1, n_construct_aop
         if (D_xy(aop, 1) .lt. min_diff) min_diff = D_xy(aop, 1)
      enddo
!
      call mem%dealloc(D_xy, n_construct_aop, 1)
!
      write(output%unit, '(/t2, a)')'- Testing the Cholesky decomposition decomposition electronic repulsion integrals:'
!
      write(output%unit, '(/t6, a71, e12.4)')'Maximal difference between approximate and actual diagonal:            ', max_diff
      write(output%unit, '(t6, a71, e12.4)')'Minimal element of difference between approximate and actual diagonal: ', min_diff
      flush(output%unit)
!
   end subroutine cholesky_vecs_diagonal_test_eri_cd_solver
!
!
   subroutine full_test_cholesky_vecs_cd_eri_solver(solver, system)
!!
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system), intent(in) :: system
!
      logical, dimension(:), allocatable :: construct_sp
!
      integer ::  A, B, C, D, w, x, y, z, wx, yz, wx_full, size_AB, AB_offset, J, I, xy
      integer ::  sp, sig_aop_counter, yx, yz_full, n_construct_sp, n_construct_aop
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
      real(dp), dimension(:,:), allocatable :: L_xy, L_xy_full, g_wxyz, g_ABCD
      integer, dimension(:,:), allocatable :: index_full
!
      character(len=40) :: line
!
!     Read significant sp info 
!
      call disk%open_file(solver%diagonal_info_construct, 'read')
!
      rewind(solver%diagonal_info_construct%unit)
!
      read(solver%diagonal_info_construct%unit) n_construct_sp, n_construct_aop
!
      allocate(construct_sp(solver%n_sp))
!
      read(solver%diagonal_info_construct%unit) construct_sp
!
      call disk%close_file(solver%diagonal_info_construct)
!
      call disk%open_file(solver%cholesky_ao_vectors, 'read')
      call disk%open_file(solver%cholesky_ao_vectors_info, 'read')
!
      rewind(solver%cholesky_ao_vectors%unit)
      rewind(solver%cholesky_ao_vectors_info%unit)
!
      call mem%alloc(L_xy, n_construct_aop, solver%n_cholesky)
      call mem%alloc(L_xy_full, solver%n_ao**2, solver%n_cholesky)
!
      L_xy = zero
      L_xy_full = zero
!
      read(solver%cholesky_ao_vectors_info%unit, '(a40)') line
!
      AB_offset = 0
!
      do while (trim(line) .ne. 'DONE')
!
!       Calculate difference between actual and approximate diagonal
!
         read(line, *) size_AB
!
         do J = 1, solver%n_cholesky
!
            read(solver%cholesky_ao_vectors%unit) (L_xy(AB_offset + I, J), I = 1, size_AB)
!
         enddo
!  
         AB_offset = AB_offset + size_AB
!
         read(solver%cholesky_ao_vectors_info%unit, *) line
!
      enddo
!
      AB_offset = 0
      sp        = 0 
      sig_aop_counter = 0
!
      call mem%alloc(index_full, n_construct_aop, 2)
!
      do B = 1, solver%n_s
!
         B_interval = system%shell_limits(B)
!
         do A = B, solver%n_s
!
            A_interval = system%shell_limits(A)
!
            sp = sp + 1
!
            if (construct_sp(sp)) then
!
              if (A .ne. B) then 
!
                 do x = A_interval%first, A_interval%last
                    do y = B_interval%first, B_interval%last
! 
                        sig_aop_counter = AB_offset + &
                              A_interval%size*(y - B_interval%first) + x - A_interval%first + 1
!
                        index_full(sig_aop_counter, 1) = x
                        index_full(sig_aop_counter, 2) = y
! 
                    enddo
                 enddo
!
              else
!
                 do x = A_interval%first, A_interval%last
                    do y = B_interval%first, B_interval%last
!
                           xy = (max(x - A_interval%first + 1,y - B_interval%first + 1)&
                              *(max(x - A_interval%first + 1, y - B_interval%first + 1)-3)/2) &
                              + x - A_interval%first + 1 + y - B_interval%first + 1
!
                           sig_aop_counter = AB_offset + xy

                           index_full(sig_aop_counter, 1) = x
                           index_full(sig_aop_counter, 2) = y
! 
                    enddo
                 enddo
!
              endif
!
               AB_offset = AB_offset + get_size_sp(A_interval, B_interval)       
!
            endif
!
         enddo
      enddo
!
      do i = 1, n_construct_aop
!
         xy = solver%n_ao*(index_full(i, 2) - 1) + index_full(i, 1)
         yx = solver%n_ao*(index_full(i, 1) - 1) + index_full(i, 2)
!
         L_xy_full(xy, 1:solver%n_cholesky) = L_xy(i, 1:solver%n_cholesky)
         L_xy_full(yx, 1:solver%n_cholesky) = L_xy(i, 1:solver%n_cholesky)
!
      enddo
!
     call mem%dealloc(L_xy, n_construct_aop, solver%n_cholesky)
!
     call mem%alloc(g_wxyz, solver%n_ao**2, solver%n_ao**2)
!
     g_wxyz = zero
!
     do A = 1, solver%n_s
        do B = 1, solver%n_s
           do C = 1, solver%n_s
              do D = 1, solver%n_s
!
                 A_interval = system%shell_limits(A)
                 B_interval = system%shell_limits(B)
                 C_interval = system%shell_limits(C)
                 D_interval = system%shell_limits(D)
!
                 call mem%alloc(g_ABCD, A_interval%size*B_interval%size, C_interval%size*D_interval%size) 
!
                 call system%ao_integrals%construct_ao_g_wxyz(g_ABCD, A, B, C, D)
!
                 do w = 1, A_interval%size
                    do x = 1, B_interval%size
                       do y = 1, C_interval%size
                          do z = 1, D_interval%size 
!
                             wx = A_interval%size*(x-1) + w
                             wx_full = solver%n_ao*(x + B_interval%first - 2) + w + A_interval%first - 1
!
                             yz = C_interval%size*(z-1) + y
                             yz_full = solver%n_ao*(z + D_interval%first - 2) + y + C_interval%first - 1
!
                             g_wxyz(wx_full, yz_full) = g_ABCD(wx, yz)
!
                          enddo
                       enddo
                    enddo
                 enddo
!
                 call mem%dealloc(g_ABCD, A_interval%size*B_interval%size, C_interval%size*D_interval%size)
!
              enddo
           enddo
        enddo
     enddo
!
     call dgemm('N','T',              &
                solver%n_ao**2,       &
                solver%n_ao**2,       &
                solver%n_cholesky,    &
                -one,                 &
                L_xy_full,            &
                solver%n_ao**2,       &
                L_xy_full,            &
                solver%n_ao**2,       &
                one,                  &           
                g_wxyz,               &
                solver%n_ao**2)
!
      write(output%unit, '(t6, a71, e12.4)')'Maximal difference between approximate and actual ERI-matrix:          ', &
                               get_abs_max(g_wxyz, solver%n_ao**4)
!
      call mem%dealloc(g_wxyz, solver%n_ao**2, solver%n_ao**2)
      call mem%dealloc(L_xy_full, solver%n_ao**2, solver%n_cholesky)
!
      call disk%close_file(solver%cholesky_ao_vectors)
      call disk%close_file(solver%cholesky_ao_vectors_info)
!
   end subroutine  full_test_cholesky_vecs_cd_eri_solver
!
!
   integer function get_size_sp(A_interval, B_interval)
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
   integer function get_sp_from_shells(s1, s2, n_s)
!!
!!    Get shell pair from shells,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    OBS THIS WILL ONLY WORK IF s1 >= s2
!!
      implicit none
!
      integer, intent(in) :: s1, s2, n_s
      integer :: A, B
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
               elseif (line(1:8) == 'batches:') then
!
                  read(line(9:100), '(i5)') solver%n_batches
!
               elseif (line(1:10) == 'qualified:') then
!
                  read(line(11:100), '(i5)') solver%max_qual
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
   subroutine construct_mo_cholesky_vecs_cd_eri_solver(solver, system, n_mo, orbital_coefficients)
!!
!!    Construct MO Cholesky vectors
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Reads AO Cholesky vectors, transforms them to MO basis and writes them to file.
!!    MO Cholesky vectors L_pq_J are stored packed (q .le. p) on direct access file with record length n_cholesky
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      type(molecular_system), intent(in) :: system
!
      integer, intent(in) :: n_mo
!
      real(dp), dimension(solver%n_ao, n_mo), intent(in) :: orbital_coefficients
!
      logical, dimension(:), allocatable :: construct_sp
!
      integer :: n_construct_sp, n_construct_aop, size_AB, A, B, I, J, AB_offset, AB_offset_full
      integer :: x, y, p, q, xy, xy_packed, sp, pq
!
      type(interval) :: A_interval, B_interval
!
      real(dp), dimension(:,:), allocatable :: L_xy, L_xy_full, L_J_pq, temp, L_pq_J, L_pq
!
      type(file) :: cholesky_mo_vectors_seq
!
      character(len=40) :: line
!
      type(batching_index) :: batch_q
!
      integer ::  current_q_batch, req0, req1, pq_rec
!
      call cholesky_mo_vectors_seq%init('cholesky_mo_temp', 'sequential', 'unformatted')
      call solver%cholesky_mo_vectors%init('cholesky_mo_vectors', 'direct', 'unformatted', dp*solver%n_cholesky)
!
      call disk%open_file(solver%cholesky_ao_vectors, 'read')
      call disk%open_file(cholesky_mo_vectors_seq, 'readwrite')
      call disk%open_file(solver%cholesky_mo_vectors, 'write')
      call disk%open_file(solver%cholesky_ao_vectors_info, 'read')
!
      rewind(solver%cholesky_ao_vectors%unit)
      rewind(solver%cholesky_ao_vectors_info%unit)
!
!     Read significant sp info 
!
      call disk%open_file(solver%diagonal_info_construct, 'read')
!
      rewind(solver%diagonal_info_construct%unit)
!
      read(solver%diagonal_info_construct%unit) n_construct_sp, n_construct_aop
!
      allocate(construct_sp(solver%n_sp))
!
      read(solver%diagonal_info_construct%unit) construct_sp
!
      call disk%close_file(solver%diagonal_info_construct)
!
      do J = 1, solver%n_cholesky
!
         rewind(solver%cholesky_ao_vectors%unit)
         rewind(solver%cholesky_ao_vectors_info%unit)
!
         read(solver%cholesky_ao_vectors_info%unit, '(a40)') line
!
         AB_offset = 0
!
         call mem%alloc(L_xy, n_construct_aop, 1)
         call mem%alloc(L_xy_full, solver%n_ao, solver%n_ao)
!
         L_xy = zero
         L_xy_full = zero 
!
         do while (trim(line) .ne. 'DONE')
!
!          Calculate difference between actual and approximate diagonal
!
            read(line, *) size_AB
!
!           Empty reads 
!
            do I = 1, J - 1
!
               read(solver%cholesky_ao_vectors%unit)
!
            enddo
!
            read(solver%cholesky_ao_vectors%unit) (L_xy(AB_offset + I, 1), I = 1, size_AB)
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

         AB_offset = 0
         AB_offset_full = 0
         sp = 0
!
         do B = 1, solver%n_s
!
            B_interval = system%shell_limits(B)
!
            do A = B, solver%n_s
!
               A_interval = system%shell_limits(A)
!
               sp = sp + 1
!
               if (construct_sp(sp)) then
!
                  if (A .ne. B) then 
!
                     do x = 1, A_interval%size
                        do y = 1, B_interval%size
!  
                            xy = A_interval%size*(y-1) + x
                            L_xy_full(x + A_interval%first - 1, y + B_interval%first - 1) = L_xy(xy + AB_offset, 1)
                            L_xy_full(y + B_interval%first - 1, x + A_interval%first - 1) = L_xy(xy + AB_offset, 1)
!  
                        enddo
                     enddo
!
                  else
!
                     do x = 1, A_interval%size
                        do y = 1, B_interval%size
!  
                            xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
                            L_xy_full(x + A_interval%first - 1, y + B_interval%first - 1) = L_xy(xy_packed + AB_offset, 1)
                            L_xy_full(y + B_interval%first - 1, x + A_interval%first - 1) = L_xy(xy_packed + AB_offset, 1)
!  
                        enddo
                     enddo
!
                  endif
!
                  AB_offset = AB_offset + get_size_sp(A_interval, B_interval)
!
               endif
!
            enddo
         enddo
!
         call mem%dealloc(L_xy, 1, n_construct_aop)
!
!       Transform the AO vectors to form the Cholesky MO vectors
!
        call mem%alloc(temp, solver%n_ao, n_mo)
!
        call dgemm('N','N',               &
                    solver%n_ao,          &
                    n_mo,                 &
                    solver%n_ao,          &
                    one,                  &
                    L_xy_full,            &
                    solver%n_ao,          &
                    orbital_coefficients, &
                    solver%n_ao,          &
                    zero,                 &
                    temp,                 &
                    solver%n_ao)
!
         call mem%alloc(L_J_pq, n_mo, n_mo)
         call mem%dealloc(L_xy_full, solver%n_ao, solver%n_ao)
!
        call dgemm('T','N',               &
                    n_mo,                 &
                    n_mo,                 &
                    solver%n_ao,          &
                    one,                  &
                    orbital_coefficients, &
                    solver%n_ao,          &
                    temp,                 &
                    solver%n_ao,          &
                    zero,                 &
                    L_J_pq,               &
                    n_mo)
!
         call mem%dealloc(temp, solver%n_ao, n_mo)
!
        write(cholesky_mo_vectors_seq%unit) L_J_pq
!
        call mem%dealloc(L_J_pq, n_mo, n_mo)
!
      enddo
!
!     Read L_pq_J in batches over q
!
      req0 = n_mo**2
      req1 = (n_mo)*(solver%n_cholesky)
!
!     Initialize batching variable
!
      call batch_q%init(n_mo)
      call mem%batch_setup(batch_q, req0, req1)
!
!     Loop over the q-batches
!
      do current_q_batch = 1, batch_q%num_batches
!
         rewind(cholesky_mo_vectors_seq%unit)
!
!        Determine limits of the q-batch
!
         call batch_q%determine_limits(current_q_batch)
!
         call mem%alloc(L_pq_J, (batch_q%length)*(n_mo), solver%n_cholesky)
!    
         call mem%alloc(L_pq, n_mo, n_mo)
!
         do J = 1, solver%n_cholesky
!
            read(cholesky_mo_vectors_seq%unit) L_pq
!
            do q = 1, batch_q%length
               do p = 1, n_mo
!
                  pq = n_mo*(q - 1) + p
!
                  L_pq_J(pq, J) = L_pq(p, q + batch_q%first - 1)
!
               enddo
            enddo
!
         enddo
!
         call mem%dealloc(L_pq, n_mo, n_mo)
!
         do p = 1, n_mo
            do q = batch_q%first, batch_q%last
!
               pq_rec = (max(p, q)*(max(p, q)-3)/2) + p + q
               pq = n_mo*(q - batch_q%first) + p
               write(solver%cholesky_mo_vectors%unit, rec=pq_rec) (L_pq_J(pq, J), J = 1, solver%n_cholesky)
!
            enddo
         enddo
!
         call mem%dealloc(L_pq_J, (batch_q%length)*(n_mo), solver%n_cholesky)
!
      enddo ! batch_q
!
      call disk%close_file(solver%cholesky_ao_vectors, 'delete')
      call disk%close_file(solver%cholesky_mo_vectors)
      call disk%close_file(cholesky_mo_vectors_seq)
      call disk%close_file(solver%cholesky_ao_vectors_info, 'delete')
!
      call disk%open_file(cholesky_mo_vectors_seq, 'read')
      call disk%close_file(cholesky_mo_vectors_seq, 'delete')
!
   end subroutine  construct_mo_cholesky_vecs_cd_eri_solver
!
!
   subroutine print_banner_eri_cd_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1)
      call long_string_print(solver%description2,'(/t3,a/)')
      call long_string_print(solver%description3)
!
   end subroutine print_banner_eri_cd_solver
!
!
   subroutine print_settings_eri_cd_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd_solver) :: solver
!
      write(output%unit, '(/t3, a)')      '- Cholesky decomposition settings:'
!
      write(output%unit, '(/t6, a21, e10.2)') 'Target threshold is: ', solver%threshold
      write(output%unit, '( t6, a21, e10.2)') 'Span factor:         ', solver%span
      write(output%unit, '( t6, a21, 4x, i6)')'Max qual:            ', solver%max_qual
      flush(output%unit)
!
   end subroutine print_settings_eri_cd_solver
!
!
end module eri_cd_solver_class
