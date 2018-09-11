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
   use reordering
!
   implicit none
!
!  Class definition
!
   type :: mo_integral_tool
!
      logical, private :: eri_file           = .false.
      logical, private :: eri_t1_file        = .false. 
      logical, private :: cholesky_file      = .true.
      logical, private :: cholesky_t1_file   = .false.
!
      integer(i15) :: n_J
!
      type(file) :: cholesky_mo
      type(file) :: cholesky_mo_t1
!
      integer(i15), private :: n_o 
      integer(i15), private :: n_v 
!
      integer(i15), private :: n_mo 
!
   contains
!
      procedure :: prepare => prepare_mo_integral_tool
  !    procedure :: cleanup => cleanup_mo_integral_tool
!
      procedure :: need_t1 => need_t1_mo_integral_tool
!
      procedure :: read_cholesky => read_cholesky_mo_integral_tool
      procedure :: read_cholesky_t1 => read_cholesky_t1_mo_integral_tool
!
      procedure :: read_cholesky_ia => read_cholesky_ia_mo_integral_tool
!
      procedure :: construct_cholesky_ij => construct_cholesky_ij_mo_integral_tool
      procedure :: construct_cholesky_ab => construct_cholesky_ab_mo_integral_tool
      procedure :: construct_cholesky_ai => construct_cholesky_ai_mo_integral_tool
!
      procedure :: read_cholesky_ai_t1 => read_cholesky_ai_t1_mo_integral_tool
      procedure :: read_cholesky_ij_t1 => read_cholesky_ij_t1_mo_integral_tool
      procedure :: read_cholesky_ab_t1 => read_cholesky_ab_t1_mo_integral_tool
!
      procedure :: set_full_index => set_full_index_mo_integral_tool
!
      procedure :: read_ovov => read_ovov_mo_integral_tool
!
      procedure :: construct_oooo => construct_oooo_mo_integral_tool
      procedure :: construct_ooov => construct_ooov_mo_integral_tool
      procedure :: construct_oovo => construct_oovo_mo_integral_tool
      procedure :: construct_ovoo => construct_ovoo_mo_integral_tool
      procedure :: construct_vooo => construct_vooo_mo_integral_tool
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
      integer(i15), intent(in) :: n_J
      integer(i15), intent(in) :: n_o
      integer(i15), intent(in) :: n_v
!
      integrals%n_J = n_J
!
      integrals%n_o  = n_o 
      integrals%n_v  = n_v
      integrals%n_mo = n_o + n_v
!
      call  integrals%cholesky_mo%init('cholesky_mo_vectors', 'direct', 'unformatted', dp*n_J)
!
   end subroutine prepare_mo_integral_tool
!
!
   logical function need_t1_mo_integral_tool(integrals)
!!
!!    Need t1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Returns true if neither t1-transformed integrals nor t1-transformed cholesky
!!    vectors are on file.
!!
      implicit none
!
      class(mo_integral_tool) :: integrals 
!
      need_t1_mo_integral_tool = ((.not. integrals%eri_t1_file) .and. (.not. integrals%cholesky_t1_file))
!
   end function need_t1_mo_integral_tool
!
!
   subroutine read_cholesky_mo_integral_tool(integrals, L_pq_J, first_p, last_p, first_q, last_q)
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
      integer(i15), intent(in) :: first_p, last_p
      integer(i15), intent(in) :: first_q, last_q
!
      real(dp), dimension((last_p - first_p + 1)*(last_q - first_q + 1), integrals%n_J) :: L_pq_J
!
      integer(i15) :: p, q, pq, pq_rec, J, dim_p, dim_q
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
!
      call disk%open_file(integrals%cholesky_mo, 'read')
!
      do p = 1, dim_p
         do q = 1, dim_q
!
            pq = dim_p*(q - 1) + p
            pq_rec = (max(p + first_p - 1, q + first_q - 1)*(max(p + first_p - 1, q + first_q - 1)-3)/2) &
                         + (p + first_p - 1) + (q + first_q - 1)
!
            read(integrals%cholesky_mo%unit, rec=pq_rec) (L_pq_J(pq, J), J = 1, integrals%n_J)
!
         enddo
      enddo

!
      call disk%close_file(integrals%cholesky_mo)
!
   end subroutine read_cholesky_mo_integral_tool
!
!
   subroutine read_cholesky_t1_mo_integral_tool(integrals, L_pq_J, first_p, last_p, first_q, last_q)
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
      integer(i15), intent(in) :: first_p, last_p
      integer(i15), intent(in) :: first_q, last_q
!
      real(dp), dimension((last_p - first_p + 1)*(last_q - first_q + 1), integrals%n_J) :: L_pq_J
!
      integer(i15) :: p, q, pq, pq_rec, J, dim_p, dim_q
!
      dim_p = last_p - first_p + 1
      dim_q = last_q - first_q + 1
!
      call disk%open_file(integrals%cholesky_mo_t1, 'read')
!
      do p = 1, dim_p
         do q = 1, dim_q
!
            pq = dim_p*(q - 1) + p
!
            read(integrals%cholesky_mo_t1%unit, rec=pq) (L_pq_J(pq, J), J = 1, integrals%n_J)
!
         enddo
      enddo

!
      call disk%close_file(integrals%cholesky_mo_t1)
!
   end subroutine read_cholesky_t1_mo_integral_tool
!
!
   subroutine read_cholesky_ia_mo_integral_tool(integrals, L_ia_J, first_i, last_i, first_a, last_a)
!!
!!    Read Cholesky ia 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!  
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_a, last_a
!
      real(dp), dimension(:,:) :: L_ia_J
!
      integer(i15) :: full_first_a, full_last_a 
      integer(i15) :: full_first_i, full_last_i
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
      call integrals%read_cholesky(L_ia_J, full_first_i, full_last_i, full_first_a, full_last_a)
!
   end subroutine read_cholesky_ia_mo_integral_tool
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
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_j, last_j
!
      real(dp), dimension(:,:) :: L_ij_J 
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: t1 
!
      real(dp), dimension(:,:), allocatable :: L_iJ_j_term
      real(dp), dimension(:,:), allocatable :: L_iJ_a
      real(dp), dimension(:,:), allocatable :: L_ia_J
!
      integer(i15) :: full_first_i, full_last_i 
      integer(i15) :: full_first_j, full_last_j 
!
      integer(i15) :: i_length, j_length          
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
      call mem%alloc(L_ia_J, i_length*(integrals%n_v), integrals%n_J)
!
!     Read the untransformed Cholesky vectors
!
      call integrals%read_cholesky(L_ij_J, full_first_i, full_last_i, full_first_j, full_last_j)
      call integrals%read_cholesky(L_ia_J, full_first_i, full_last_i, integrals%n_o + 1, integrals%n_mo)
!
!     Compute and add t1-transformed term, L_iJ_j sum_a t_aj L_ia_J 
!
      call mem%alloc(L_iJ_a, i_length*(integrals%n_J), integrals%n_v)
      call sort_123_to_132(L_ia_J, L_iJ_a, i_length, integrals%n_v, integrals%n_J)
      call mem%dealloc(L_ia_J, i_length*(integrals%n_v), integrals%n_J)
!
      call mem%alloc(L_iJ_j_term, i_length*(integrals%n_J), j_length)
!
      call dgemm('N','N',                   &
                  i_length*(integrals%n_J), &
                  j_length,                 &
                  integrals%n_v,            &
                  one,                      &
                  L_iJ_a,                   &
                  i_length*(integrals%n_J), &
                  t1(1, full_first_j),      &
                  integrals%n_v,            &
                  zero,                     &
                  L_iJ_j_term,              &
                  i_length*(integrals%n_J))
!
     call add_132_to_123(one, L_iJ_j_term, L_ij_J, i_length, j_length, integrals%n_J)
!
     call mem%dealloc(L_iJ_j_term, i_length*(integrals%n_J), j_length)
     call mem%dealloc(L_iJ_a, i_length*(integrals%n_J), integrals%n_v)
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
      real(dp), dimension(:,:) :: L_ab_J
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: t1 
!
      integer(i15), intent(in) :: first_a, last_a
      integer(i15), intent(in) :: first_b, last_b
!
      integer(i15) :: full_first_a, full_last_a 
      integer(i15) :: full_first_b, full_last_b
!
      integer(i15) :: red_first_a 
!
      real(dp), dimension(:,:), allocatable :: L_ib_J
      real(dp), dimension(:,:), allocatable :: L_Jb_i
      real(dp), dimension(:,:), allocatable :: L_Jb_a
!
      integer(i15) :: b_length, a_length
!
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
      call integrals%set_full_index(full_first_b, 'f', 'v', first_b)
!
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
      call integrals%set_full_index(full_last_b, 'l', 'v', last_b)      
!
      red_first_a = full_first_a - integrals%n_o
!
      a_length = full_last_a - full_first_a + 1
      b_length = full_last_b - full_first_b + 1
!
!     Set L_ib_J = L_ib^J and L_ab_J^T1 = L_ab_J 
!
      call mem%alloc(L_ib_J, (integrals%n_o)*b_length, integrals%n_J)
      call integrals%read_cholesky(L_ib_J, 1, integrals%n_o, full_first_b, full_last_b)
!
      call integrals%read_cholesky(L_ab_J, full_first_a, full_last_a, full_first_b, full_last_b)
!
!     Reorder L_ib^J as L_Jb_i 
!
      call mem%alloc(L_Jb_i, (integrals%n_J)*b_length, integrals%n_o)
      call sort_123_to_321(L_ib_J, L_Jb_i, integrals%n_o, b_length, integrals%n_J)
      call mem%dealloc(L_ib_J, (integrals%n_o)*b_length, integrals%n_J)
!
!     Calculate and add t1-transformed term, - sum_i t_ai L_ib_J
!
      call mem%alloc(L_Jb_a, (integrals%n_J)*b_length, a_length)
!
      call dgemm('N','T',                   &
                  (integrals%n_J)*b_length, &
                  a_length,                 &
                  integrals%n_o,            &
                  -one,                     &
                  L_Jb_i,                   &
                  (integrals%n_J)*b_length, &
                  t1(red_first_a, 1),       &
                  integrals%n_v,            &
                  zero,                     &
                  L_Jb_a,                   &
                  b_length*(integrals%n_J))
!
      call add_321_to_123(one, L_Jb_a, L_ab_J, a_length, b_length, integrals%n_J)
!
      call mem%dealloc(L_Jb_a, (integrals%n_J)*b_length, a_length)
      call mem%dealloc(L_Jb_i, (integrals%n_J)*b_length, integrals%n_o)      
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
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_a, last_a
!
      real(dp), dimension(:, :) :: L_ai_J
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
      real(dp), dimension(:, :), allocatable :: L_bj_J, L_ji_J, L_ba_J, X_j_iJ, X_i_jJ, X_i_aJ
!
      integer(i15) :: full_first_a, full_last_a, length_a 
      integer(i15) :: full_first_i, full_last_i, length_i
!
      integer(i15) :: i, J, a, ai, aJ, current_a_batch, required
!
      type(batching_index) :: batch_a
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
      length_i = full_last_i - full_first_i + 1
      length_a = full_last_a - full_first_a + 1

      call integrals%read_cholesky(L_ai_J, full_first_a, full_last_a, full_first_i, full_last_i)
!
      call mem%alloc(L_bj_J, (integrals%n_v)*(integrals%n_o), 1)
      call integrals%read_cholesky(L_bj_J, 1, (integrals%n_v), 1, (integrals%n_o))
!
      call mem%alloc(X_i_jJ, length_i, (integrals%n_o)*(integrals%n_J))
!
      call dgemm('T', 'N',                                  &
                  length_i,                                 &
                  (integrals%n_o)*(integrals%n_J),          &
                  integrals%n_v,                            &
                  one,                                      &
                  t1(1, full_first_i),                      &   ! t_b_i
                  integrals%n_v,                            &
                  L_bj_J,                                   &   ! L_b_jJ
                  integrals%n_v,                            &
                  zero,                                     &
                  X_i_jJ,                                   &
                  length_i)
!
      call mem%dealloc(L_bj_J, (integrals%n_v)*(integrals%n_o), 1)
!
      call mem%alloc(X_j_iJ,(integrals%n_o), length_i*(integrals%n_J))
!
      call sort_123_to_213(X_i_jJ, X_j_iJ, length_i, (integrals%n_o), (integrals%n_J))
!
      call mem%dealloc(X_i_jJ, length_i, (integrals%n_o)*(integrals%n_J))
!
      call dgemm('N', 'N',                                  &
                  length_a,                                 &
                  (length_i)*(integrals%n_J),               &
                  integrals%n_o,                            &
                  -one,                                     &
                  t1(full_first_a, 1),                      & ! t_a_j
                  integrals%n_v,                            &
                  X_j_iJ,                                   &
                  integrals%n_o,                            &
                  one,                                      &
                  L_ai_J,                                   & ! L_a_iJ
                  length_a)
!
      call mem%dealloc(X_j_iJ, (integrals%n_o), length_i*(integrals%n_J))
!
      call mem%alloc(L_ji_J, (integrals%n_o)*length_i, (integrals%n_J))
      call integrals%read_cholesky(L_ji_J, 1, (integrals%n_o), full_first_i, full_last_i)
!
      call dgemm('N', 'N',                                  &
                  length_a,                                 &
                  (length_i)*(integrals%n_J),               &
                  integrals%n_o,                            &
                  -one,                                     &
                  t1(full_first_a, 1),                      & ! t_a_j
                  integrals%n_v,                            &
                  L_ji_J,                                   & ! L_j_iJ
                  integrals%n_o,                            &
                  one,                                      &
                  L_ai_J,                                   & ! L_a_iJ
                  length_a)     
!
      call mem%dealloc(L_ji_J, (integrals%n_o)*length_i, (integrals%n_J))
!
      call batch_a%init(length_a)
!
      required = length_a*(integrals%n_v)*(integrals%n_J) &
               + length_i*(length_a)*(integrals%n_J)
!
      call mem%num_batch(batch_a, required)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch) 
!
         call mem%alloc(L_ba_J, batch_a%length*(integrals%n_v), (integrals%n_J)) 
!
         call integrals%read_cholesky(L_ba_J, 1, (integrals%n_v), &
                        batch_a%first + (integrals%n_o), batch_a%last + (integrals%n_o))
!
         call mem%alloc(X_i_aJ, length_i, (batch_a%length)*(integrals%n_J)) 
!
         call dgemm('T', 'N',                               &
                  length_i,                                 &
                  (batch_a%length)*(integrals%n_J),         &
                  integrals%n_v,                            &
                  one,                                      &
                  t1,                                       & ! t_b_i
                  integrals%n_v,                            &
                  L_ba_J,                                   & ! L_b_aJ
                  integrals%n_v,                            &
                  one,                                      &
                  X_i_aJ,                                   & 
                  length_i)  
!        
         call mem%dealloc(L_ba_J, batch_a%length*(integrals%n_v), (integrals%n_J))  
!
!$omp parallel do &
!$omp private(i, a, J, ai, aJ) &
!$omp shared(L_ai_J, X_i_aJ, batch_a)
         do i = 1, length_i
            do a = 1, batch_a%length
               do J = 1, (integrals%n_J)
!
                  ai = length_a*(i - 1) + a + batch_a%first - 1
                  aJ = batch_a%length*(J - 1) + a
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + X_i_aJ(i, aJ)
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(X_i_aJ, length_i, (batch_a%length)*(integrals%n_J)) 
!
       enddo
!
   end subroutine construct_cholesky_ai_mo_integral_tool
!
!
   subroutine read_cholesky_ai_t1_mo_integral_tool(integrals, L_ai_J, first_a, last_a, first_i, last_i)
!!
!!    Read Cholesky ai T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Read T1-transformed mo cholesky L_ai^J vectors from file 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!  
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_a, last_a
!
      real(dp), dimension(:, :) :: L_ai_J
!
      integer(i15) :: full_first_a, full_last_a, length_a 
      integer(i15) :: full_first_i, full_last_i, length_i
!
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
      call integrals%read_cholesky_t1(L_ai_J, full_first_a, full_last_a, full_first_i, full_last_i)
!
   end subroutine read_cholesky_ai_t1_mo_integral_tool
!
!
   subroutine read_cholesky_ij_t1_mo_integral_tool(integrals, L_ij_J, first_i, last_i, first_j, last_j)
!!
!!    Read Cholesky ij T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Read T1-transformed mo cholesky L_ij^J vectors from file 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!  
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_j, last_j
!
      real(dp), dimension(:, :) :: L_ij_J
!
      integer(i15) :: full_first_j, full_last_j 
      integer(i15) :: full_first_i, full_last_i
!
!
      call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
      call integrals%set_full_index(full_first_j, 'f', 'v', first_j)
!
      call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
      call integrals%set_full_index(full_last_j, 'l', 'v', last_j)
!
      call integrals%read_cholesky_t1(L_ij_J, full_first_i, full_last_i, full_first_j, full_last_j)
!
   end subroutine read_cholesky_ij_t1_mo_integral_tool
!
!
   subroutine read_cholesky_ab_t1_mo_integral_tool(integrals, L_ab_J, first_a, last_a, first_b, last_b)
!!
!!    Read Cholesky ab T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Read T1-transformed mo cholesky L_ab^J vectors from file 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!  
      integer(i15), intent(in) :: first_a, last_a
      integer(i15), intent(in) :: first_b, last_b
!
      real(dp), dimension(:, :) :: L_ab_J
!
      integer(i15) :: full_first_b, full_last_b 
      integer(i15) :: full_first_a, full_last_a
!
      call integrals%set_full_index(full_first_a, 'f', 'o', first_a)
      call integrals%set_full_index(full_first_b, 'f', 'v', first_b)
!
      call integrals%set_full_index(full_last_a, 'l', 'o', last_a)
      call integrals%set_full_index(full_last_b, 'l', 'v', last_b)
!
      call integrals%read_cholesky_t1(L_ab_J, full_first_a, full_last_a, full_first_b, full_last_b)
!
   end subroutine read_cholesky_ab_t1_mo_integral_tool
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
      integer(i15), intent(inout) :: ind 
!
      character(len=1), intent(in) :: orb_space 
      character(len=1), intent(in) :: pos  
!
      integer(i15), optional, intent(in) :: red_ind 
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
   subroutine construct_oooo_mo_integral_tool(integrals, g_ijkl, t1, first_i, last_i, first_j, last_j, &
                                                            first_k, last_k, first_l, last_l, index_restrictions)
!!
!!    Construct oooo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_ijkl
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_k, last_k
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_l, last_l
!
      integer(i15) :: length_i, length_k, length_j, length_l
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ij_J 
      real(dp), dimension(:,:), allocatable :: L_kl_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_k = last_k - first_k + 1
      length_l = last_l - first_l + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_ijkl from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_ijkl from Cholesky vectors
!
         if (index_restrictions) then ! dim_ij ≠ dim_kl in general
!
            call mem%alloc(L_ij_J, length_i*length_j, integrals%n_J)
            call mem%alloc(L_kl_J, length_k*length_l, integrals%n_J)
!
            call integrals%construct_cholesky_ij(L_ij_J, t1, first_i, last_i, first_j, last_j)
            call integrals%construct_cholesky_ij(L_kl_J, t1, first_k, last_k, first_l, last_l)
!
            call dgemm('N', 'T',           &
                        length_i*length_j, &
                        length_k*length_l, &
                        integrals%n_J,     &
                        one,               &
                        L_ij_J,            &
                        length_i*length_j, &
                        L_kl_J,            &
                        length_k*length_l, &
                        zero,              &
                        g_ijkl)
!
            call mem%dealloc(L_ij_J, length_i*length_j, integrals%n_J)
            call mem%dealloc(L_kl_J, length_k*length_l, integrals%n_J)
!
         else ! dim_ij = dim_kl
!
            call mem%alloc(L_ij_J, length_i*length_j, integrals%n_J)
!
            call integrals%construct_cholesky_ij(L_ij_J, t1, first_i, last_i, first_j, last_j)
!
            call dgemm('N', 'T',           &
                        length_i*length_j, &
                        length_i*length_j, &
                        integrals%n_J,     &
                        one,               &
                        L_ij_J,            &
                        length_i*length_j, &
                        L_ij_J,            &
                        length_i*length_j, & 
                        zero,              &
                        g_ijkl)
!
            call mem%dealloc(L_ij_J, length_i*length_j, integrals%n_J)
!
         endif 
!
      endif 
!
   end subroutine construct_oooo_mo_integral_tool
!
!
   subroutine construct_ooov_mo_integral_tool(integrals, g_ijka, t1, first_i, last_i, first_j, last_j, &
                                                            first_k, last_k, first_a, last_a, index_restrictions)
!!
!!    Construct ooov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_ijka
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_k, last_k
      integer(i15), optional, intent(in) :: first_j, last_j
      integer(i15), optional, intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_k, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ij_J 
      real(dp), dimension(:,:), allocatable :: L_ka_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_k = last_k - first_k + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_ijka from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_ijka
!
         call mem%alloc(L_ij_J, length_i*length_j, integrals%n_J)
         call mem%alloc(L_ka_J, length_k*length_a, integrals%n_J)
!
         call integrals%construct_cholesky_ij(L_ij_J, t1, first_i, last_i, first_j, last_j)
         call integrals%read_cholesky_ia(L_ka_J, first_k, last_k, first_a, last_a)
!
         call dgemm('N', 'T',           &
                     length_i*length_j, &
                     length_k*length_a, &
                     integrals%n_J,     &
                     one,               &
                     L_ij_J,            &
                     length_i*length_j, &
                     L_ka_J,            &
                     length_k*length_a, &
                     zero,              &
                     g_ijka)
!
         call mem%dealloc(L_ij_J, length_i*length_j, integrals%n_J)
         call mem%dealloc(L_ka_J, length_k*length_a, integrals%n_J)
!
      endif 
!
   end subroutine construct_ooov_mo_integral_tool
!
!
   subroutine construct_oovo_mo_integral_tool(integrals, g_ijak, t1, first_i, last_i, first_j, last_j, &
                                                            first_a, last_a, first_k, last_k, index_restrictions)
!!
!!    Construct oovo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_ijak
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_k, last_k
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_k, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ij_J 
      real(dp), dimension(:,:), allocatable :: L_ak_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_k = last_k - first_k + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_ijak from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_ijak
!
         call mem%alloc(L_ij_J, length_i*length_j, integrals%n_J)
         call mem%alloc(L_ak_J, length_k*length_a, integrals%n_J)
!
         call integrals%construct_cholesky_ij(L_ij_J, t1, first_i, last_i, first_j, last_j)
         call integrals%construct_cholesky_ai(L_ak_J, t1, first_a, last_a, first_k, last_k)
!
         call dgemm('N', 'T',           &
                     length_i*length_j, &
                     length_k*length_a, &
                     integrals%n_J,     &
                     one,               &
                     L_ij_J,            &
                     length_i*length_j, &
                     L_ak_J,            &
                     length_k*length_a, &
                     zero,              &
                     g_ijak)
!
         call mem%dealloc(L_ij_J, length_i*length_j, integrals%n_J)
         call mem%dealloc(L_ak_J, length_k*length_a, integrals%n_J)
!
      endif 
!
   end subroutine construct_oovo_mo_integral_tool
!
!
   subroutine construct_ovoo_mo_integral_tool(integrals, g_iajk, t1, first_i, last_i, first_a, last_a, &
                                                            first_j, last_j, first_k, last_k, index_restrictions)
!!
!!    Construct ovoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_iajk
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_k, last_k
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_k, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_jk_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_k = last_k - first_k + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_iajk from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_iajk
!
         call mem%alloc(L_ia_J, length_i*length_a, integrals%n_J)
         call mem%alloc(L_jk_J, length_k*length_j, integrals%n_J)
!
         call integrals%construct_cholesky_ij(L_jk_J, t1, first_j, last_j, first_k, last_k)
         call integrals%read_cholesky_ia(L_ia_J, first_i, last_i, first_a, last_a)
!
         call dgemm('N', 'T',           &
                     length_i*length_a, &
                     length_k*length_j, &
                     integrals%n_J,     &
                     one,               &
                     L_ia_J,            &
                     length_i*length_a, &
                     L_jk_J,            &
                     length_k*length_j, &
                     zero,              &
                     g_iajk)
!
         call mem%dealloc(L_ia_J, length_i*length_a, integrals%n_J)
         call mem%dealloc(L_jk_J, length_k*length_j, integrals%n_J)
!
      endif 
!
   end subroutine construct_ovoo_mo_integral_tool
!
!
   subroutine construct_vooo_mo_integral_tool(integrals, g_aijk, t1, first_a, last_a, first_i, last_i, &
                                                            first_j, last_j, first_k, last_k, index_restrictions)
!!
!!    Construct ovoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_aijk
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_k, last_k
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_k, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_jk_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_k = last_k - first_k + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_aijk from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_aijk
!
         call mem%alloc(L_ai_J, length_i*length_a, integrals%n_J)
         call mem%alloc(L_jk_J, length_k*length_j, integrals%n_J)
!
         call integrals%construct_cholesky_ij(L_jk_J, t1, first_j, last_j, first_k, last_k)
         call integrals%construct_cholesky_ai(L_ai_J, t1, first_a, last_a, first_i, last_i)
!
         call dgemm('N', 'T',           &
                     length_i*length_a, &
                     length_k*length_j, &
                     integrals%n_J,     &
                     one,               &
                     L_ai_J,            &
                     length_a*length_i, &
                     L_jk_J,            &
                     length_j*length_k, &
                     zero,              &
                     g_aijk)
!
         call mem%dealloc(L_ai_J, length_i*length_a, integrals%n_J)
         call mem%dealloc(L_jk_J, length_k*length_j, integrals%n_J)
!
      endif 
!
   end subroutine construct_vooo_mo_integral_tool
!
!
   subroutine construct_vvoo_mo_integral_tool(integrals, g_abij, t1, first_a, last_a, first_b, last_b, &
                                                            first_i, last_i, first_j, last_j, index_restrictions)
!!
!!    Construct vvoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_abij
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_b, last_b
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_b, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ab_J 
      real(dp), dimension(:,:), allocatable :: L_ij_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_b = last_b - first_b + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_abij from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_abij
!
         call mem%alloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%alloc(L_ij_J, length_i*length_j, integrals%n_J)
!
         call integrals%construct_cholesky_ij(L_ij_J, t1, first_i, last_i, first_j, last_j)
         call integrals%construct_cholesky_ab(L_ab_J, t1, first_a, last_a, first_b, last_b)
!
         call dgemm('N', 'T',           &
                     length_a*length_b, &
                     length_i*length_j, &
                     integrals%n_J,     &
                     one,               &
                     L_ab_J,            &
                     length_a*length_b, &
                     L_ij_J,            &
                     length_i*length_j, &
                     zero,              &
                     g_abij)
!
         call mem%dealloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%dealloc(L_ij_J, length_i*length_j, integrals%n_J)
!
      endif 
!
   end subroutine construct_vvoo_mo_integral_tool
!
!
!
   subroutine construct_oovv_mo_integral_tool(integrals, g_ijab, t1, first_i, last_i, first_j, last_j, &
                                                            first_a, last_a, first_b, last_b, index_restrictions)
!!
!!    Construct oovv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_ijab
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_b, last_b
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_b, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ab_J 
      real(dp), dimension(:,:), allocatable :: L_ij_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_b = last_b - first_b + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_ijab from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_ijab
!
         call mem%alloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%alloc(L_ij_J, length_i*length_j, integrals%n_J)
!
         call integrals%construct_cholesky_ij(L_ij_J, t1, first_i, last_i, first_j, last_j)
         call integrals%construct_cholesky_ab(L_ab_J, t1, first_a, last_a, first_b, last_b)
!
         call dgemm('N', 'T',           &
                     length_i*length_j, &
                     length_a*length_b, &
                     integrals%n_J,     &
                     one,               &
                     L_ij_J,            &
                     length_i*length_j, &
                     L_ab_J,            &
                     length_a*length_b, &
                     zero,              &
                     g_ijab)
!
         call mem%dealloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%dealloc(L_ij_J, length_i*length_j, integrals%n_J)
!
      endif 
!
   end subroutine construct_oovv_mo_integral_tool
!
!
   subroutine construct_voov_mo_integral_tool(integrals, g_aijb, t1, first_a, last_a, first_i, last_i, &
                                                            first_j, last_j, first_b, last_b, index_restrictions)
!!
!!    Construct voov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_aijb
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_b, last_b
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_b, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_jb_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_b = last_b - first_b + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_aijb from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_aijb
!
         call mem%alloc(L_ai_J, length_a*length_i, integrals%n_J)
         call mem%alloc(L_jb_J, length_b*length_j, integrals%n_J)
!
         call integrals%construct_cholesky_ai(L_ai_J, t1, first_a, last_a, first_i, last_i)
         call integrals%read_cholesky_ia(L_jb_J, first_j, last_j, first_b, last_b)
!
         call dgemm('N', 'T',           &
                     length_a*length_i, &
                     length_b*length_j, &
                     integrals%n_J,     &
                     one,               &
                     L_ai_J,            &
                     length_i*length_a, &
                     L_jb_J,            &
                     length_j*length_b, &
                     zero,              &
                     g_aijb)
!
         call mem%dealloc(L_ai_J, length_a*length_i, integrals%n_J)
         call mem%dealloc(L_jb_J, length_b*length_j, integrals%n_J)
!
      endif 
!
   end subroutine construct_voov_mo_integral_tool
!
!
   subroutine construct_ovvo_mo_integral_tool(integrals, g_iabj, t1, first_i, last_i, first_a, last_a, &
                                                            first_b, last_b, first_j, last_j, index_restrictions)
!!
!!    Construct voov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_iabj
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_b, last_b
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_a, last_a
!
      integer(i15) :: length_i, length_b, length_j, length_a
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_bj_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_b = last_b - first_b + 1
      length_a = last_a - first_a + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_iabj from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_iabj
!
         call mem%alloc(L_ia_J, length_a*length_i, integrals%n_J)
         call mem%alloc(L_bj_J, length_b*length_j, integrals%n_J)
!
         call integrals%construct_cholesky_ai(L_bj_J, t1, first_b, last_b, first_j, last_j)
         call integrals%read_cholesky_ia(L_ia_J, first_i, last_i, first_a, last_a)
!
         call dgemm('N', 'T',           &
                     length_a*length_i, &
                     length_b*length_j, &
                     integrals%n_J,     &
                     one,               &
                     L_ia_J,            &
                     length_i*length_a, &
                     L_bj_J,            &
                     length_j*length_b, &
                     zero,              &
                     g_iabj)
!
         call mem%dealloc(L_ia_J, length_a*length_i, integrals%n_J)
         call mem%dealloc(L_bj_J, length_b*length_j, integrals%n_J)
!
      endif 
!
   end subroutine construct_ovvo_mo_integral_tool
!
!
   subroutine read_ovov_mo_integral_tool(integrals, g_iajb, first_i, last_i, first_a, last_a, &
                                                            first_j, last_j, first_b, last_b, index_restrictions)
!!
!!    Read ovov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Although this is named read, note that it is used also when t1 transformed 
!!    integrals are needed, since there is in this case - iajb - no contribution from the
!!    t1 amplitudes: the MO integrals are equal to the t1-transformed MO integrals.
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_iajb
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_a, last_a
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_b, last_b
!
      integer(i15) :: length_i, length_a, length_j, length_b
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_jb_J 
!
      length_i = last_i - first_i + 1
      length_a = last_a - first_a + 1
      length_j = last_j - first_j + 1
      length_b = last_b - first_b + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_iajb from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_iajb from Cholesky vectors (occ-vir, so not necessary to consider t1)
!
         if (index_restrictions) then ! dim_ia ≠ dim_jb in general
!
            call mem%alloc(L_ia_J, length_i*length_a, integrals%n_J)
            call mem%alloc(L_jb_J, length_j*length_b, integrals%n_J)
!
            call integrals%read_cholesky_ia(L_ia_J, first_i, last_i, first_a, last_a)
            call integrals%read_cholesky_ia(L_jb_J, first_j, last_j, first_b, last_b)
!
            call dgemm('N', 'T',           &
                        length_i*length_a, &
                        length_j*length_b, &
                        integrals%n_J,     &
                        one,               &
                        L_ia_J,            &
                        length_i*length_a, &
                        L_jb_J,            &
                        length_j*length_b, &
                        zero,              &
                        g_iajb)
!
            call mem%dealloc(L_ia_J, length_i*length_a, integrals%n_J)
            call mem%dealloc(L_jb_J, length_j*length_b, integrals%n_J)
!
         else ! dim_ia = dim_jb 
!
            call mem%alloc(L_ia_J, length_i*length_a, integrals%n_J)
!
            call integrals%read_cholesky_ia(L_ia_J, first_i, last_i, first_a, last_a)
!
            call dgemm('N', 'T',           &
                        length_i*length_a, &
                        length_i*length_a, &
                        integrals%n_J,     &
                        one,               &
                        L_ia_J,            &
                        length_i*length_a, &
                        L_ia_J,            &
                        length_i*length_a, & 
                        zero,              &
                        g_iajb)
!
            call mem%dealloc(L_ia_J, length_i*length_a, integrals%n_J)
!
         endif 
!
      endif 
!
   end subroutine read_ovov_mo_integral_tool
!
!
   subroutine construct_vovo_mo_integral_tool(integrals, g_aibj, t1, first_a, last_a, first_i, last_i, &
                                                            first_b, last_b, first_j, last_j, index_restrictions)
!!
!!    Construct vovo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_aibj
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), intent(in) :: first_i, last_i
      integer(i15), intent(in) :: first_a, last_a
      integer(i15), intent(in) :: first_j, last_j
      integer(i15), intent(in) :: first_b, last_b
!
      integer(i15) :: length_i, length_a, length_j, length_b
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_bj_J 
!
      length_i = last_i - first_i + 1
      length_j = last_j - first_j + 1
      length_a = last_a - first_a + 1
      length_b = last_b - first_b + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_aibj from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_aibj from Cholesky vectors
!
         if (index_restrictions) then ! dim_ia ≠ dim_jb in general
!
            call mem%alloc(L_ai_J, length_a*length_i, integrals%n_J)
            call mem%alloc(L_bj_J, length_b*length_j, integrals%n_J)
!
            call integrals%construct_cholesky_ai(L_ai_J, t1, first_a, last_a, first_i, last_i)
            call integrals%construct_cholesky_ai(L_bj_J, t1, first_b, last_b, first_j, last_j)
!
            call dgemm('N', 'T',           &
                        length_a*length_i, &
                        length_b*length_j, &
                        integrals%n_J,     &
                        one,               &
                        L_ai_J,            &
                        length_i*length_a, &
                        L_bj_J,            &
                        length_b*length_j, &
                        zero,              &
                        g_aibj)
!
            call mem%dealloc(L_ai_J, length_a*length_i, integrals%n_J)
            call mem%dealloc(L_bj_J, length_b*length_j, integrals%n_J)
!
         else ! dim_ai = dim_bj
!
            call mem%alloc(L_ai_J, length_a*length_i, integrals%n_J)
!
            call integrals%construct_cholesky_ai(L_ai_J, t1, first_a, last_a, first_i, last_i)
!
            call dgemm('N', 'T',           &
                        length_i*length_a, &
                        length_i*length_a, &
                        integrals%n_J,     &
                        one,               &
                        L_ai_J,            &
                        length_i*length_a, &
                        L_ai_J,            &
                        length_i*length_a, & 
                        zero,              &
                        g_aibj)
!
            call mem%dealloc(L_ai_J, length_i*length_a, integrals%n_J)
!
         endif 
!
      endif 
!
   end subroutine construct_vovo_mo_integral_tool
!
!
   subroutine construct_vvvo_mo_integral_tool(integrals, g_abci, t1, first_a, last_a, first_b, last_b, &
                                                            first_c, last_c, first_i, last_i, index_restrictions)
!!
!!    Construct vvvo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_abci
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_a, last_a
      integer(i15), optional, intent(in) :: first_b, last_b
      integer(i15), optional, intent(in) :: first_c, last_c
!
      integer(i15) :: length_i, length_a, length_b, length_c
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ab_J 
      real(dp), dimension(:,:), allocatable :: L_ci_J 
!
      length_i = last_i - first_i + 1
      length_a = last_a - first_a + 1
      length_b = last_b - first_b + 1
      length_c = last_c - first_c + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_abci from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_abci
!
         call mem%alloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%alloc(L_ci_J, length_c*length_i, integrals%n_J)
!
         call integrals%construct_cholesky_ab(L_ab_J, t1, first_a, last_a, first_b, last_b)
         call integrals%construct_cholesky_ai(L_ci_J, t1, first_c, last_c, first_i, last_i)
!
         call dgemm('N', 'T',           &
                     length_a*length_b, &
                     length_c*length_i, &
                     integrals%n_J,     &
                     one,               &
                     L_ab_J,            &
                     length_a*length_b, &
                     L_ci_J,            &
                     length_c*length_i, &
                     zero,              &
                     g_abci)
!
         call mem%dealloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%dealloc(L_ci_J, length_c*length_i, integrals%n_J)
!
      endif 
!
   end subroutine construct_vvvo_mo_integral_tool
!
!
   subroutine construct_vvov_mo_integral_tool(integrals, g_abic, t1, first_a, last_a, first_b, last_b, &
                                                            first_i, last_i, first_c, last_c, index_restrictions)
!!
!!    Construct vvov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_abic
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_a, last_a
      integer(i15), optional, intent(in) :: first_b, last_b
      integer(i15), optional, intent(in) :: first_c, last_c
!
      integer(i15) :: length_i, length_a, length_b, length_c
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ab_J 
      real(dp), dimension(:,:), allocatable :: L_ic_J 
!
      length_i = last_i - first_i + 1
      length_a = last_a - first_a + 1
      length_b = last_b - first_b + 1
      length_c = last_c - first_c + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_abic from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_abic
!
         call mem%alloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%alloc(L_ic_J, length_c*length_i, integrals%n_J)
!
         call integrals%construct_cholesky_ab(L_ab_J, t1, first_a, last_a, first_b, last_b)
         call integrals%read_cholesky_ia(L_ic_J, first_i, last_i, first_c, last_c)
!
         call dgemm('N', 'T',           &
                     length_a*length_b, &
                     length_c*length_i, &
                     integrals%n_J,     &
                     one,               &
                     L_ab_J,            &
                     length_a*length_b, &
                     L_ic_J,            &
                     length_c*length_i, &
                     zero,              &
                     g_abic)
!
         call mem%dealloc(L_ab_J, length_a*length_b, integrals%n_J)
         call mem%dealloc(L_ic_J, length_c*length_i, integrals%n_J)
!
      endif 
!
   end subroutine construct_vvov_mo_integral_tool
!
!
   subroutine construct_vovv_mo_integral_tool(integrals, g_aibc, t1, first_a, last_a, first_i, last_i, &
                                                            first_b, last_b, first_c, last_c, index_restrictions)
!!
!!    Construct vovv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_aibc
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_a, last_a
      integer(i15), optional, intent(in) :: first_b, last_b
      integer(i15), optional, intent(in) :: first_c, last_c
!
      integer(i15) :: length_i, length_a, length_b, length_c
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ai_J 
      real(dp), dimension(:,:), allocatable :: L_bc_J 
!
      length_i = last_i - first_i + 1
      length_a = last_a - first_a + 1
      length_b = last_b - first_b + 1
      length_c = last_c - first_c + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_aibc from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_aibc
!
         call mem%alloc(L_ai_J, length_a*length_i, integrals%n_J)
         call mem%alloc(L_bc_J, length_c*length_b, integrals%n_J)
!
         call integrals%construct_cholesky_ai(L_ai_J, t1, first_a, last_a, first_i, last_i)
         call integrals%construct_cholesky_ab(L_bc_J, t1, first_b, last_b, first_c, last_c)
!
         call dgemm('N', 'T',           &
                     length_a*length_i, &
                     length_c*length_b, &
                     integrals%n_J,     &
                     one,               &
                     L_ai_J,            &
                     length_a*length_i, &
                     L_bc_J,            &
                     length_c*length_b, &
                     zero,              &
                     g_aibc)
!
         call mem%dealloc(L_ai_J, length_a*length_i, integrals%n_J)
         call mem%dealloc(L_bc_J, length_c*length_b, integrals%n_J)
!
      endif 
!
   end subroutine construct_vovv_mo_integral_tool
!
!
   subroutine construct_ovvv_mo_integral_tool(integrals, g_iabc, t1, first_i, last_i, first_a, last_a, &
                                                            first_b, last_b, first_c, last_c, index_restrictions)
!!
!!    Construct ovvv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:), intent(inout) :: g_iabc
!
      real(dp), dimension(integrals%n_v, integrals%n_o) :: t1
!
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_a, last_a
      integer(i15), optional, intent(in) :: first_b, last_b
      integer(i15), optional, intent(in) :: first_c, last_c
!
      integer(i15) :: length_i, length_a, length_b, length_c
!
      logical, intent(in) :: index_restrictions
!
      real(dp), dimension(:,:), allocatable :: L_ia_J 
      real(dp), dimension(:,:), allocatable :: L_bc_J 
!
      length_i = last_i - first_i + 1
      length_a = last_a - first_a + 1
      length_b = last_b - first_b + 1
      length_c = last_c - first_c + 1
!
      if (integrals%eri_file .and. .not. index_restrictions) then 
!
!        Coming soon: read full g_aibc from file
!
         call output%error_msg('reading full eri integrals from file not yet supported!')
!
      else
!
!        Construct g_aibc
!
         call mem%alloc(L_ia_J, length_a*length_i, integrals%n_J)
         call mem%alloc(L_bc_J, length_c*length_b, integrals%n_J)
!
         call integrals%read_cholesky_ia(L_ia_J, first_a, last_a, first_i, last_i)
         call integrals%construct_cholesky_ab(L_bc_J, t1, first_b, last_b, first_c, last_c)
!
         call dgemm('N', 'T',           &
                     length_a*length_i, &
                     length_c*length_b, &
                     integrals%n_J,     &
                     one,               &
                     L_ia_J,            &
                     length_a*length_i, &
                     L_bc_J,            &
                     length_c*length_b, &
                     zero,              &
                     g_iabc)
!
         call mem%dealloc(L_ia_J, length_a*length_i, integrals%n_J)
         call mem%dealloc(L_bc_J, length_c*length_b, integrals%n_J)
!
      endif 
!
   end subroutine construct_ovvv_mo_integral_tool
!
!
end module mo_integral_tool_class

