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
      integer(i15) :: n_cholesky
!
      type(file) :: cholesky_mo
      type(file) :: cholesky_mo_t1
!
   contains
!
  !    procedure :: prepare => prepare_mo_integral_tool
  !    procedure :: cleanup => cleanup_mo_integral_tool
!
      procedure :: need_t1 => need_t1_mo_integral_tool
!
      procedure :: read_cholesky => read_cholesky_mo_integral_tool
!
   end type mo_integral_tool
!
!
contains
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
   subroutine read_cholesky_mo_integral_tool(integrals, L_pq_J, first_p, dim_p, first_q, dim_q)
!!
!!    Read mo cholesky vectors
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Reads cholesky vectors L_pq_J for mo indices 
!!    [first_p, first_p + dim_p - 1] and [first_q, first_q + dim_q - 1]
!!
      implicit none
!
      class(mo_integral_tool) :: integrals
!
      integer(i15), intent(in) :: first_p, dim_p
      integer(i15), intent(in) :: first_q, dim_q
!
      real(dp), dimension(dim_p*dim_q, integrals%n_cholesky) :: L_pq_J
!
      integer(i15) :: p, q, pq, pq_rec, J
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
            read(integrals%cholesky_mo%unit, rec=pq_rec) (L_pq_J(pq, J), J = 1, integrals%n_cholesky)
!
         enddo
      enddo

!
      call disk%close_file(integrals%cholesky_mo)
!
   end subroutine read_cholesky_mo_integral_tool
!
!
end module mo_integral_tool_class

