submodule (ccs_class) integrals
!
!
!
contains
!
!
   module subroutine get_oo_oo_ccs(wf, integral_type, x_oo_oo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
!
!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_oo_oo
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_oo_electronic_repulsion(x_oo_oo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_oo_oo'
!
      endif
!
   end subroutine get_oo_oo_ccs
!
!
   module subroutine get_vo_vo_ccs(wf, integral_type, x_vo_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vo_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vo_vo_electronic_repulsion(x_vo_vo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_vo_vo'
!
      endif
!
   end subroutine get_vo_vo_ccs
!
!
   module subroutine get_ov_vo_ccs(wf, integral_type, x_ov_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_ov_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_vo_electronic_repulsion(x_ov_vo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_ov_vo'
!
      endif
!
   end subroutine get_ov_vo_ccs
!
!
   module subroutine get_vo_ov_ccs(wf, integral_type, x_vo_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vo_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_vo_electronic_repulsion(x_vo_ov,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_vo_ov'
!
      endif
!
   end subroutine get_vo_ov_ccs
!
!
   module subroutine get_ov_vv_ccs(wf, integral_type, x_ov_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_ov_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_vv_electronic_repulsion(x_ov_vv,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_ov_vv'
!
      endif
!
   end subroutine get_ov_vv_ccs
!
!
   module subroutine get_vv_ov_ccs(wf, integral_type, x_vv_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vv_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_ov_electronic_repulsion(x_vv_ov,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_vv_ov'
!
      endif
!
   end subroutine get_vv_ov_ccs
!
!
   module subroutine get_vo_vv_ccs(wf, integral_type, x_vo_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vo_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vo_vv_electronic_repulsion(x_vo_vv,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_vo_vv'
!
      endif
!
   end subroutine get_vo_vv_ccs
!
!
   module subroutine get_vv_vo_ccs(wf, integral_type, x_vv_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vv_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_vo_electronic_repulsion(x_vv_vo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_vv_vo'
!
      endif
!
   end subroutine get_vv_vo_ccs
!
!
   module subroutine get_vv_vv_ccs(wf, integral_type, x_vv_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vv_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_vo_electronic_repulsion(x_vv_vv,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_vv_vv'
!
      endif
!
   end subroutine get_vv_vv_ccs
!
!
   module subroutine get_vo_vo_electronic_repulsion_ccs(wf, x_vo_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vo_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_bj_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ai_J, length_1*length_2, wf%n_J)
      call allocator(L_bj_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ai(L_bj_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ai_J,              &
                  length_1*length_2,   &
                  L_bj_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_vo_vo,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ai_J, length_1*length_2, wf%n_J)
      call deallocator(L_bj_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Get T1-transformed Cholesky vector
!
      call wf%get_cholesky_ai(L_ai_J)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_J,              &
                  one,                 &
                  L_ai_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  L_ai_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  x_vo_vo,             &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vo_vo_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_vo_vo_electronic_repulsion_ccs
!
!
module subroutine get_ov_vo_electronic_repulsion_ccs(wf, x_ov_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_ov_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_bj_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ia_J, length_1*length_2, wf%n_J)
      call allocator(L_bj_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ia(L_ia_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ai(L_bj_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ia_J,              &
                  length_1*length_2,   &
                  L_bj_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_ov_vo,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ia_J, length_1*length_2, wf%n_J)
      call deallocator(L_bj_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ia_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call allocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ia(L_ia_J)
      call wf%get_cholesky_ai(L_bj_J)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_J,              &
                  one,                 &
                  L_ia_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  L_bj_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  x_ov_vo,             &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ia_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call deallocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_ov_vo_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_ov_vo_electronic_repulsion_ccs
!
!
module subroutine get_vo_ov_electronic_repulsion_ccs(wf, x_vo_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vo_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_jb_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ai_J, length_1*length_2, wf%n_J)
      call allocator(L_jb_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ia(L_jb_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ai_J,              &
                  length_1*length_2,   &
                  L_jb_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_vo_ov,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ai_J, length_1*length_2, wf%n_J)
      call deallocator(L_jb_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call allocator(L_jb_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ai(L_ai_J)
      call wf%get_cholesky_ia(L_jb_J)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_J,              &
                  one,                 &
                  L_ai_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  L_jb_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  x_vo_ov,             &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call deallocator(L_jb_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vo_ov_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_vo_ov_electronic_repulsion_ccs
!
!
   module subroutine get_ov_vv_electronic_repulsion_ccs(wf, x_ov_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_ov_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_bc_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ia_J, length_1*length_2, wf%n_J)
      call allocator(L_bc_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ia(L_ia_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ab(L_bc_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ia_J,              &
                  length_1*length_2,   &
                  L_bc_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_ov_vv,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ia_J, length_1*length_2, wf%n_J)
      call deallocator(L_bc_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ia_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call allocator(L_bc_J, (wf%n_v)**2, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ia(L_ia_J)
      call wf%get_cholesky_ab(L_bc_J, 1, wf%n_v, 1, wf%n_v)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)**2,         &
                  wf%n_J,              &
                  one,                 &
                  L_ia_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  L_bc_J,              &
                  (wf%n_v)**2,         &
                  zero,                &
                  x_ov_vv,             &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ia_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call deallocator(L_bc_J, (wf%n_v)**2, wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_ov_vv_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_ov_vv_electronic_repulsion_ccs
!
!
   module subroutine get_vv_ov_electronic_repulsion_ccs(wf, x_vv_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vv_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_ic_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ab_J, length_1*length_2, wf%n_J)
      call allocator(L_ic_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ia(L_ic_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  length_1*length_2,   &
                  L_ic_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_vv_ov,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ab_J, length_1*length_2, wf%n_J)
      call deallocator(L_ic_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
      call allocator(L_ic_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, 1, wf%n_v)
      call wf%get_cholesky_ia(L_ic_J)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)**2,         &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  (wf%n_v)**2,         &
                  L_ic_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  x_vv_ov,             &
                  (wf%n_v)**2)
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
      call deallocator(L_ic_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vv_ov_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_vv_ov_electronic_repulsion_ccs
!
!
   module subroutine get_vo_vv_electronic_repulsion_ccs(wf, x_vo_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vo_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_bc_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ai_J, length_1*length_2, wf%n_J)
      call allocator(L_bc_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ab(L_bc_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ai_J,              &
                  length_1*length_2,   &
                  L_bc_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_vo_vv,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ai_J, length_1*length_2, wf%n_J)
      call deallocator(L_bc_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call allocator(L_bc_J, (wf%n_v)**2, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ai(L_ai_J)
      call wf%get_cholesky_ia(L_bc_J, 1, wf%n_v, 1, wf%n_v)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)**2,         &
                  wf%n_J,              &
                  one,                 &
                  L_ai_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  L_bc_J,              &
                  (wf%n_v)**2,         &
                  zero,                &
                  x_vo_vv,             &
                  (wf%n_v)*(wf%n_o))
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)      
      call deallocator(L_bc_J, (wf%n_v)**2, wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vo_vv_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_vo_vv_electronic_repulsion_ccs
!
!
   module subroutine get_vv_vo_electronic_repulsion_ccs(wf, x_vv_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vv_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_ci_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ab_J, length_1*length_2, wf%n_J)
      call allocator(L_ci_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ai(L_ci_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  length_1*length_2,   &
                  L_ci_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_vv_vo,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ab_J, length_1*length_2, wf%n_J)
      call deallocator(L_ci_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
      call allocator(L_ci_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, 1, wf%n_v)
      call wf%get_cholesky_ia(L_ci_J)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)**2,         &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  (wf%n_v)**2,         &
                  L_ci_J,              &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  x_vv_vo,             &
                  (wf%n_v)**2)
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
      call deallocator(L_ci_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vv_vo_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_vv_vo_electronic_repulsion_ccs
!
!
!
   module subroutine get_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vv_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_cd_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
   if (present(index1_first) .and. present(index1_first)      &
      .and. present(index2_first) .and. present(index2_first) &
      .and. present(index3_first) .and. present(index3_first) &
      .and. present(index4_first) .and. present(index4_first) ) then
!
!     Optional arguments are pressent, we are either batching or running MLCC calculation
!
!     Sanity check here!!
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Alllocate Cholesky vectors
!
      call allocator(L_ab_J, length_1*length_2, wf%n_J)
      call allocator(L_cd_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ab(L_cd_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  length_1*length_2,   &
                  L_cd_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_vv_vv,             &
                  length_1*length_2)
!
!     Deallocate Cholesky vectors
!
      call deallocator(L_ab_J, length_1*length_2, wf%n_J)
      call deallocator(L_cd_J, length_3*length_4, wf%n_J)
!
   elseif ( .not. (present(index1_first) .and. present(index1_first) &
         .and. present(index2_first) .and. present(index2_first)     &
         .and. present(index3_first) .and. present(index3_first)     &
         .and. present(index4_first) .and. present(index4_first) )) then
!
!     No optional arguments passed
!
!     Alllocate Cholesky vector
!
      call allocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, 1, wf%n_v)
!
!     Construct integral
!
      call dgemm('N', 'T',             &
                  (wf%n_v)**2,         &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  (wf%n_v)**2,         &
                  L_ab_J,              &
                  (wf%n_v)**2,         &
                  zero,                &
                  x_vv_vv,             &
                  (wf%n_v)**2)
!
!     Deallocate Cholesky vector
!
      call deallocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
!
   else
!
!     Something wrong in subroutine call
!
      write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vv_vv_electronic_repulsion'
      stop
!
   endif
!
   end subroutine get_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine get_oo_oo_electronic_repulsion_ccs(wf, x_oo_oo,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_oo_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_kl_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      if (present(index1_first) .and. & 
          present(index1_last)  .and. &
          present(index2_first) .and. &
          present(index2_last)  .and. &
          present(index3_first) .and. &
          present(index3_last)  .and. & 
          present(index4_first) .and. &
          present(index4_last)) then
!
!        Optional arguments are pressent, we are either batching or running MLCC calculation
!
!        Sanity check here!!
!
!        Lengths
!
         length_1 = index1_last - index1_first + 1
         length_2 = index2_last - index2_first + 1
         length_3 = index3_last - index3_first + 1
         length_4 = index4_last - index4_first + 1
!
!        Alllocate Cholesky vectors
!
         call allocator(L_ij_J, length_1*length_2, wf%n_J)
         call allocator(L_kl_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ij(L_kl_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
!
         call dgemm('N', 'T',           &
                     length_1*length_2, &
                     length_3*length_4, &
                     wf%n_J,            &
                     one,               &
                     L_ij_J,            &
                     length_1*length_2, &
                     L_kl_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_oo_oo,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ij_J, length_1*length_2, wf%n_J)
         call deallocator(L_kl_J, length_3*length_4, wf%n_J)
!
      elseif ( .not. (present(index1_first) .and. & 
                      present(index1_first) .and. &
                      present(index2_first) .and. &
                      present(index2_first) .and. &
                      present(index3_first) .and. &
                      present(index3_first) .and. &
                      present(index4_first) .and. &
                      present(index4_first)) ) then
!
!        No optional arguments passed
!
!        Alllocate Cholesky vector
!
         call allocator(L_ij_J, (wf%n_o)**2, wf%n_J)      
!
!        Get entire T1-transformed Cholesky vector
!
         call wf%get_cholesky_ij(L_ij_J)
!
!        Construct integral
!
         call dgemm('N', 'T',  &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  wf%n_J,      &
                  one,         &
                  L_ij_J,      &
                  (wf%n_o)**2, &
                  L_ij_J,      &
                  (wf%n_o)**2, &
                  zero,        &
                  x_oo_oo,     &
                  (wf%n_o)**2)
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J)      
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vv_vv_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_oo_oo_electronic_repulsion_ccs
!
!
end submodule integrals