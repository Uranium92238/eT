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
!
      procedure :: get_cholesky_ia => get_cholesky_ia_mo_integral_tool
      procedure :: get_cholesky_ij => get_cholesky_ij_mo_integral_tool
!
      procedure :: set_full_index => set_full_index_mo_integral_tool
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
   subroutine get_cholesky_ia_mo_integral_tool(integrals, L_ia_J, first_i, last_i, first_a, last_a)
!!
!!    Get Cholesky ia 
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
   end subroutine get_cholesky_ia_mo_integral_tool
!
!
   subroutine get_cholesky_ij_mo_integral_tool(integrals, L_ij_J, t1, first_i, last_i, first_j, last_j)
!!
!!    Get Cholesky ij
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
   end subroutine get_cholesky_ij_mo_integral_tool
!
!
   subroutine get_cholesky_ab_mo_integral_tool(integrals, L_ab_J, t1, first_a, last_a, first_b, last_b)
!!
!!    Get Cholesky ab 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Computes 
!!
!!       L_ab_J_T1 = L_ab_J - sum_i t_ai L_ib_J
!! 
!!    and saves the result in L_ab_J.
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!
      real(dp), dimension(:,:) :: L_ab_J
!
      real(dp), dimension(integrals%n_v, integrals%n_o), intent(in) :: t1 
!
      integer(i15), optional, intent(in) :: first_a, last_a
      integer(i15), optional, intent(in) :: first_b, last_b
!
      integer(i15) :: full_first_a, full_last_a 
      integer(i15) :: full_first_b, full_last_b
!
      call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
      call integrals%set_full_index(full_first_b, 'f', 'v', first_b)
!
      call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
      call integrals%set_full_index(full_last_b, 'l', 'v', last_b)      
!
      
!
   end subroutine get_cholesky_ab_mo_integral_tool
!
!
   subroutine get_cholesky_ai_mo_integral_tool(integrals, L_ai_J, first_a, last_a, first_i, last_i)
!!
!!    Get Cholesky ai
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Not done!!!
!!
      implicit none 
!
      class(mo_integral_tool), intent(in) :: integrals 
!  
      integer(i15), optional, intent(in) :: first_i, last_i
      integer(i15), optional, intent(in) :: first_a, last_a
!
      real(dp), dimension(:, :) :: L_ai_J
!
      integer(i15) :: full_first_a, full_last_a 
      integer(i15) :: full_first_i, full_last_i
!
    !  call integrals%set_full_index(full_first_i, 'f', 'o', first_i)
    !  call integrals%set_full_index(full_first_a, 'f', 'v', first_a)
!
    !  call integrals%set_full_index(full_last_i, 'l', 'o', last_i)
    !  call integrals%set_full_index(full_last_a, 'l', 'v', last_a)
!
    !  call integrals%read_cholesky(L_ia_J, full_first_i, full_last_i, full_first_a, full_last_a)
!
   end subroutine get_cholesky_ai_mo_integral_tool
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
end module mo_integral_tool_class

