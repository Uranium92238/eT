submodule (ccs_class) integrals
!
!!
!!    Integrals submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
!!
!!    Contains procedures for construction of general integrals, and 
!!    spesifically:
!!
!!    - electronic repulsion integrals
!!    
!!    o - occupied index
!!    v - virtual index
!
contains
!
!
   module subroutine get_oo_oo_ccs(wf, integral_type, x_oo_oo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_oo,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_oo,oo (ordered as x_oo_oo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_oo_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
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
   module subroutine get_oo_ov_ccs(wf, integral_type, x_oo_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_oo,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_oo,ov (ordered as x_oo_ov)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_oo_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_ov_electronic_repulsion(x_oo_ov,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_oo_ov'
!
      endif
!
   end subroutine get_oo_ov_ccs
!
!
   module subroutine get_ov_oo_ccs(wf, integral_type, x_ov_oo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_ov,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_ov,oo (ordered as x_ov_oo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_ov_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_oo_electronic_repulsion(x_ov_oo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*)'WARNING: unknown integral type requested from get_ov_oo'
!
      endif
!
   end subroutine get_ov_oo_ccs
!
!
   module subroutine get_oo_vo_ccs(wf, integral_type, x_oo_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_oo,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_oo,vo (ordered as x_oo_vo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_oo_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_vo_electronic_repulsion(x_oo_vo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*) 'WARNING: unknown integral type requested from get_oo_vo'
!
      endif
!
   end subroutine get_oo_vo_ccs
!
!
   module subroutine get_vo_oo_ccs(wf, integral_type, x_vo_oo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_vo,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vo,oo (ordered as x_vo_oo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vo_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vo_oo_electronic_repulsion(x_vo_oo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*) 'WARNING: unknown integral type requested from get_vo_oo'
!
      endif
!
   end subroutine get_vo_oo_ccs
!
!
   module subroutine get_oo_vv_ccs(wf, integral_type, x_oo_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_oo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_oo,vv (ordered as x_oo_vv)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_oo_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_vv_electronic_repulsion(x_oo_vv,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*) 'WARNING: unknown integral type requested from get_oo_vv'
!
      endif
!
   end subroutine get_oo_vv_ccs
!
!
   module subroutine get_vv_oo_ccs(wf, integral_type, x_vv_oo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_vv,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vv,oo (ordered as x_vv_oo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_vv_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_oo_electronic_repulsion(x_vv_oo,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*) 'WARNING: unknown integral type requested from get_vv_oo'
!
      endif
!
   end subroutine get_vv_oo_ccs
!
!
   module subroutine get_ov_ov_ccs(wf, integral_type, x_ov_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_ov,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_ov,ov (ordered as x_ov_ov)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=40) :: integral_type 
!
      real(dp), dimension(:, :) :: x_ov_ov 
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_ov_electronic_repulsion(x_ov_ov,          & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
!
         write(unit_output,*) 'WARNING: unknown integral type requested from get_ov_ov'
!
      endif
!
   end subroutine get_ov_ov_ccs
!
!
   module subroutine get_vo_vo_ccs(wf, integral_type, x_vo_vo,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Get x_vo,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vo,vo (ordered as x_vo_vo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
!!    Get x_ov,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_ov,vo (ordered as x_ov_vo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
!!    Get x_vo,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vo,ov (ordered as x_vo_ov)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
         call wf%get_vo_ov_electronic_repulsion(x_vo_ov,          & 
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
!!    Get x_ov,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_ov,vv (ordered as x_ov_vv)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
!!    Get x_vv,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vv,ov (ordered as x_vv_ov)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
!!    Get x_vo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vo,vv (ordered as x_vo_vv)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
!!    Get x_vv,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vv,vo (ordered as x_vv_vo)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
!!    Get x_vv,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct genereal two-electron integral x_vv,vv (ordered as x_vv_vv)
!!    of type integral_type.
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
         call wf%get_vv_vv_electronic_repulsion(x_vv_vv,          & 
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
!!    Get g_vo,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,vo (ordered as g_vo_vo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_ov,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_ov,vo (ordered as g_ov_vo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_vo,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,ov (ordered as g_vo_ov)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_vo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,vv (ordered as g_vo_vv)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_vv,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,ov (ordered as g_vv_ov)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_vo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,vv (ordered as g_vo_vv)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_vv,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,vo (ordered as g_vv_vo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
   if (     present(index1_first) .and. present(index1_last) &
      .and. present(index2_first) .and. present(index2_last) &
      .and. present(index3_first) .and. present(index3_last) &
      .and. present(index4_first) .and. present(index4_last) ) then
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
   elseif ( .not. (present(index1_first) .and. present(index1_last) &
             .and. present(index2_first) .and. present(index2_last) &
             .and. present(index3_first) .and. present(index3_last) &
             .and. present(index4_first) .and. present(index4_last) )) then
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
!!    Get g_vv,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,vv (ordered as g_vv_vv)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
      logical :: vvvv_on_file = .false.
!
      inquire(file='g_abcd',exist=vvvv_on_file)
      if (vvvv_on_file) then
!
!     Read x_vv_vv
!
!
!     T1-transform x_vv_vv
!
      !call wf%t1_transform_vv_vv(x_vv_vv
      else
         if (     present(index1_first) .and. present(index1_last) &
            .and. present(index2_first) .and. present(index2_last) &
            .and. present(index3_first) .and. present(index3_last) &
            .and. present(index4_first) .and. present(index4_last) ) then
!
!           Optional arguments are pressent, we are either batching or running MLCC calculation
!
!           Sanity check here!!
!
!           Lengths
!
            length_1 = index1_last - index1_first + 1
            length_2 = index2_last - index2_first + 1
            length_3 = index3_last - index3_first + 1
            length_4 = index4_last - index4_first + 1
!
!           Alllocate Cholesky vectors
!
            call allocator(L_ab_J, length_1*length_2, wf%n_J)
            call allocator(L_cd_J, length_3*length_4, wf%n_J)
!
!           Get T1-transformed Cholesky vectors
!
            call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
            call wf%get_cholesky_ab(L_cd_J, index3_first, index3_last, index4_first, index4_last)
!
!           Construct integral
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
!           Deallocate Cholesky vectors
!
            call deallocator(L_ab_J, length_1*length_2, wf%n_J)
            call deallocator(L_cd_J, length_3*length_4, wf%n_J)
!
         elseif ( .not. (present(index1_first) .and. present(index1_last) &
                   .and. present(index2_first) .and. present(index2_last) &
                   .and. present(index3_first) .and. present(index3_last) &
                   .and. present(index4_first) .and. present(index4_last) )) then
!
!           No optional arguments passed
!
!           Alllocate Cholesky vector
!
            call allocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
!
!           Get T1-transformed Cholesky vectors
!
            call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, 1, wf%n_v)
!
!           Construct integral
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
!           Deallocate Cholesky vector
!
            call deallocator(L_ab_J, (wf%n_v)**2, wf%n_J)      
!
         else
!
!           Something wrong in subroutine call
!
            write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vv_vv_electronic_repulsion'
            stop
!
         endif
      endif
!
   end subroutine get_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine read_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:) :: x_vv_vv
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last 
!
!     Temporary vector for holding parts of integral 
!
      real(dp), dimension(:,:), allocatable :: x_v ! g_a_bcd, a = 1, n_v, for given bcd
!
!     File handling integers 
!
      integer(i15) :: unit_g_abcd = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1     ! Error integer for file handling
      integer(i15) :: rec_number = -1  ! The record where g_abcd is positioned 
!
!     Indices 
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0, ab = 0, cd = 0, bcd = 0, acd = 0 
!
!     Open file containing the g_abcd integrals, ordered as
!     g_a_bcd, where the compound cd index is packed.
!
!     The compound index (b, cd_packed) determines the record number,
!     where the record includes the integrals g(a,bcd), a = 1, n_v
!
      call generate_unit_identifier(unit_g_abcd)
      open(unit=unit_g_abcd, file='g_abcd', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_v), iostat=ioerror)
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,*) 'Error: could not open file g_abcd in read_vv_vv_electronic_repulsion_ccs'
         stop
!
      endif
!
!     Allocate the integral x_v the integrals g_a_bcd for a given bcd, a = 1, n_v 
!
      call allocator(x_v, wf%n_v, 1)
      x_v = zero 
!
      if ((index1_first .ne. 1) .or. (index1_last .ne. wf%n_v)) then ! Batching over first index, a 
!
         if ((index2_first .ne. 1) .or. (index2_last .ne. wf%n_v)) then ! Batching over second index, b, as well
!
!           No simple tricks available. Read all integrals g_a_bcd, a = 1, n_v, into x_v, then through away!
!
            do d = index4_first, index4_last
               do c = index3_first, index3_last
                  do b = index2_first, index2_last
!
                     bcd = index_two(b, index_packed(c,d), wf%n_v) ! Record number 
!
!                    Read g_a_bcd, a = 1, n_v, into x_v 
!
                     read(unit_g_abcd, rec=bcd) (x_v(a, 1), a = 1, wf%n_v)
!
!                    Place the integral into x_vv_vv = g_ab_cd 
!
                     cd = index_two(c, d, index3_last - index3_first)
!
                     do a = index1_first, index1_last
!
                        ab = index_two(a, b, index1_last - index1_first)
!
                        x_vv_vv(ab, cd) = x_v(a, 1)
!
                     enddo
!
                  enddo
               enddo
            enddo      
!
         else ! Batching over first but not second index => pretend a is b to avoid reading more than necessary
!
!           Pretend first index is second index (switch a and b, such that batching is over b):
!           Read g_b_acd, b = 1, n_v, into x_v
!
            do d = index4_first, index4_last
               do c = index3_first, index3_last
                  do a = index1_first, index1_last
!
                     acd = index_two(a, index_packed(c,d), wf%n_v) ! Record number 
!
!                    Read g_b_acd, b = 1, n_v, into x_v 
!
                     read(unit_g_abcd, rec=acd) (x_v(b, 1), b = 1, wf%n_v)
!
!                    Place the integral into x_vv_vv = g_ab_cd 
!
                     cd = index_two(c, d, index3_last - index3_first)
!
                     do b = index2_first, index2_last
!
                        ab = index_two(a, b, index1_last - index1_first)
!
                        x_vv_vv(ab, cd) = x_v(b, 1) ! g_abcd 
!
                     enddo
!
                  enddo
               enddo
            enddo
!
         endif 
!
      endif 
!
!     Close file containing the g_abcd integrals 
!
      close(unit_g_abcd)
!
   end subroutine read_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine get_oo_oo_electronic_repulsion_ccs(wf, x_oo_oo,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_oo,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_oo,oo (ordered as g_oo_oo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
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
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_oo_oo_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_oo_oo_electronic_repulsion_ccs
!
!
   module subroutine get_oo_ov_electronic_repulsion_ccs(wf, x_oo_ov,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_oo,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_oo,ov (ordered as g_oo_ov)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_oo_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_ka_J
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
         call allocator(L_ka_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ia(L_ka_J, index3_first, index3_last, index4_first, index4_last)
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
                     L_ka_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_oo_ov,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ij_J, length_1*length_2, wf%n_J)
         call deallocator(L_ka_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ij_J, (wf%n_o)**2, wf%n_J) 
         call allocator(L_ka_J, (wf%n_o)*(wf%n_v), wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J)
         call wf%get_cholesky_ia(L_ka_J)
!
!        Construct integral
!
         call dgemm('N', 'T',        &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ij_J,            &
                  (wf%n_o)**2,       &
                  L_ka_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  x_oo_ov,           &
                  (wf%n_o)**2)
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J) 
         call deallocator(L_ka_J, (wf%n_o)*(wf%n_v), wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_oo_ov_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_oo_ov_electronic_repulsion_ccs
!
!
   module subroutine get_ov_oo_electronic_repulsion_ccs(wf, x_ov_oo,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_ov,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_ov,oo (ordered as g_ov_oo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_ov_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_jk_J
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
         call allocator(L_ia_J, length_1*length_2, wf%n_J)
         call allocator(L_jk_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ia(L_ia_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ij(L_jk_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
!
         call dgemm('N', 'T',           &
                     length_1*length_2, &
                     length_3*length_4, &
                     wf%n_J,            &
                     one,               &
                     L_ia_J,            &
                     length_1*length_2, &
                     L_jk_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_ov_oo,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ia_J, length_1*length_2, wf%n_J)
         call deallocator(L_jk_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J) 
         call allocator(L_jk_J, (wf%n_o)**2, wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ia(L_ia_J)
         call wf%get_cholesky_ij(L_jk_J)
!
!        Construct integral
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)**2,       &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_jk_J,            &
                  (wf%n_o)**2,       &
                  zero,              &
                  x_ov_oo,           &
                  (wf%n_o)*(wf%n_v))
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J) 
         call deallocator(L_jk_J, (wf%n_o)**2, wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_ov_oo_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_ov_oo_electronic_repulsion_ccs
!
!
   module subroutine get_oo_vo_electronic_repulsion_ccs(wf, x_oo_vo,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_oo,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_oo,vo (ordered as g_oo_vo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_oo_vo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_ak_J
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
         call allocator(L_ak_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ai(L_ak_J, index3_first, index3_last, index4_first, index4_last)
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
                     L_ak_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_oo_vo,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ij_J, length_1*length_2, wf%n_J)
         call deallocator(L_ak_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ij_J, (wf%n_o)**2, wf%n_J) 
         call allocator(L_ak_J, (wf%n_v)*(wf%n_o), wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J)
         call wf%get_cholesky_ai(L_ak_J)
!
!        Construct integral
!
         call dgemm('N', 'T',        &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ij_J,            &
                  (wf%n_o)**2,       &
                  L_ak_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  x_oo_vo,           &
                  (wf%n_o)**2)
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J) 
         call deallocator(L_ak_J, (wf%n_o)*(wf%n_v), wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_ov_oo_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_oo_vo_electronic_repulsion_ccs
!
!
   module subroutine get_vo_oo_electronic_repulsion_ccs(wf, x_vo_oo,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_vo,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vooo (ordered as g_vo_oo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vo_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_jk_J
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
         call allocator(L_ai_J, length_1*length_2, wf%n_J)
         call allocator(L_jk_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ij(L_jk_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
!
         call dgemm('N', 'T',           &
                     length_1*length_2, &
                     length_3*length_4, &
                     wf%n_J,            &
                     one,               &
                     L_ai_J,            &
                     length_1*length_2, &
                     L_jk_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_vo_oo,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ai_J, length_1*length_2, wf%n_J)
         call deallocator(L_jk_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J) 
         call allocator(L_jk_J, (wf%n_o)**2, wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ai(L_ai_J)
         call wf%get_cholesky_ij(L_jk_J)
!
!        Construct integral
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)**2,       &
                  wf%n_J,            &
                  one,               &
                  L_ai_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_jk_J,            &
                  (wf%n_o)**2,       &
                  zero,              &
                  x_vo_oo,           &
                  (wf%n_o)*(wf%n_v))
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J) 
         call deallocator(L_jk_J, (wf%n_o)**2, wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vo_oo_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_vo_oo_electronic_repulsion_ccs
!
!
   module subroutine get_oo_vv_electronic_repulsion_ccs(wf, x_oo_vv,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_oo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_oo,vv (ordered as g_oo_vv)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_oo_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_ab_J
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
         call allocator(L_ab_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ab(L_ab_J, index3_first, index3_last, index4_first, index4_last)
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
                     L_ab_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_oo_vv,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ij_J, length_1*length_2, wf%n_J)
         call deallocator(L_ab_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ij_J, (wf%n_o)**2, wf%n_J) 
         call allocator(L_ab_J, (wf%n_v)**2, wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ij(L_ij_J)
         call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, 1, wf%n_v)
!
!        Construct integral
!
         call dgemm('N', 'T',  &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  wf%n_J,      &
                  one,         &
                  L_ij_J,      &
                  (wf%n_o)**2, &
                  L_ab_J,      &
                  (wf%n_v)**2, &
                  zero,        &
                  x_oo_vv,     &
                  (wf%n_o)**2)
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J) 
         call deallocator(L_ab_J, (wf%n_v)**2, wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_oo_vv_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_oo_vv_electronic_repulsion_ccs
!
!
   module subroutine get_vv_oo_electronic_repulsion_ccs(wf, x_vv_oo,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_vv,oo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,oo (ordered as g_vv_oo)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vv_oo
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_ij_J
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
         call allocator(L_ab_J, length_1*length_2, wf%n_J)
         call allocator(L_ij_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ij(L_ij_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
!
         call dgemm('N', 'T',           &
                     length_1*length_2, &
                     length_3*length_4, &
                     wf%n_J,            &
                     one,               &
                     L_ab_J,            &
                     length_1*length_2, &
                     L_ij_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_vv_oo,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ab_J, length_1*length_2, wf%n_J)
         call deallocator(L_ij_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ab_J, (wf%n_v)**2, wf%n_J) 
         call allocator(L_ij_J, (wf%n_o)**2, wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, 1, wf%n_v)
         call wf%get_cholesky_ij(L_ij_J)
!
!        Construct integral
!
         call dgemm('N', 'T',  &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  wf%n_J,      &
                  one,         &
                  L_ab_J,      &
                  (wf%n_v)**2, &
                  L_ij_J,      &
                  (wf%n_o)**2, &
                  zero,        &
                  x_vv_oo,     &
                  (wf%n_v)**2)
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ab_J, (wf%n_v)**2, wf%n_J) 
         call deallocator(L_ij_J, (wf%n_o)**2, wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_vv_oo_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_vv_oo_electronic_repulsion_ccs
!
!
   module subroutine get_ov_ov_electronic_repulsion_ccs(wf, x_ov_ov,                & 
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!    Get g_ov,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_ov,ov (ordered as g_ov_ov)
!!
!!    Optional parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_ov_ov
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_jb_J
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
         call allocator(L_ia_J, length_1*length_2, wf%n_J)
         call allocator(L_jb_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ia(L_ia_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ia(L_jb_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
!
         call dgemm('N', 'T',           &
                     length_1*length_2, &
                     length_3*length_4, &
                     wf%n_J,            &
                     one,               &
                     L_ia_J,            &
                     length_1*length_2, &
                     L_jb_J,            &
                     length_3*length_4, &
                     zero,              &
                     x_ov_ov,           &
                     length_1*length_2)
!
!        Deallocate Cholesky vectors
!
         call deallocator(L_ia_J, length_1*length_2, wf%n_J)
         call deallocator(L_jb_J, length_3*length_4, wf%n_J)
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
!        Alllocate Cholesky vectors
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)     
!
!        Get entire T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ia(L_ia_J)
!
!        Construct integral
!
         call dgemm('N', 'T',        &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            & ! L_jb_J 
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  x_ov_ov,           &
                  (wf%n_o)*(wf%n_v))
!
!        Deallocate Cholesky vector
!
         call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)     
!
      else
!
!        Something wrong in subroutine call
!
         write(unit_output,*) 'WARNING: Some, but not all optional arguments were passed to get_ov_ov_electronic_repulsion'
         stop
!
      endif
!
   end subroutine get_ov_ov_electronic_repulsion_ccs
!
!
   module subroutine t1_transform_vv_vv_ccs(wf, g_vv_vv,                & 
                                             index1_first, index1_last, &
                                             index2_first, index2_last, &
                                             index3_first, index3_last, &
                                             index4_first, index4_last)
!!
!!       T1 transformation of g_vv_vv integrals (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Oct 2017. 
!!
!!       g_ab_cd_T1 = g_ab_cd + sum_(J) sum_(i) t_a_i * L_ib_J * L_cd_J
!!                            + sum_(J) sum_(k) t_c_k * L_kd_J * L_ab_J
!!                            + sum_(J) sum_(ki) t_a_i * t_c_k * L_kd_J * L_ib_J
!! 
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: g_vv_vv
!
      integer(i15), optional :: index1_first, index1_last
      integer(i15), optional :: index2_first, index2_last
      integer(i15), optional :: index3_first, index3_last
      integer(i15), optional :: index4_first, index4_last
!
      integer(i15) :: c = 0, d = 0, cd = 0, dc = 0
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      real(dp), dimension(:,:), allocatable :: L_ib_J, L_dk_J, L_cd_J, L_ab_J
      real(dp), dimension(:,:), allocatable :: x_ib_dk, x_ib_dc, x_ib_cd, x_ab_dk
      real(dp), dimension(:,:), allocatable :: g_ab_dc
!
      length_1 = index1_last - index1_first + 1 ! a 
      length_2 = index2_last - index2_first + 1 ! b 
      length_3 = index3_last - index3_first + 1 ! c
      length_4 = index4_last - index4_first + 1 ! d
!
!     ::Term 1::
!     - sum_(J) sum_(i) t_a_i * L_ib_J * L_cd_J 
!
      call allocator(L_ib_J, wf%n_o*length_2, wf%n_J)
      call allocator(L_cd_J, length_3*length_4, wf%n_J)
!
      call wf%read_cholesky_ia(L_ib_J, 1, wf%n_o, index2_first, index2_last)
      call wf%read_cholesky_ab(L_cd_J, index3_first, index3_last, index4_first, index4_last)
!
      call allocator(x_ib_cd, wf%n_o*length_2, length_3*length_4)
!
!     x_ib_cd = sum_(J) L_ib_J * L_cd_J
!
      call dgemm('N', 'T',             &
                  wf%n_o*length_2,     &
                  length_3*length_4,   &
                  wf%n_J,              &
                  one,                 &
                  L_ib_J,              &
                  wf%n_o*length_2,     &
                  L_cd_J,              &
                  length_3*length_4,   &
                  zero,                &
                  x_ib_cd,             &
                  wf%n_o*length_2)
!
      call deallocator(L_cd_J, length_3*length_4, wf%n_J)
!
!     g_vv_vv = g_ab_cd += sum_(i)t_a_i* x_ib_cd
!
      call dgemm('N', 'N',                      &
                  length_1,                     &
                  length_2*length_4*length_3,   &
                  wf%n_o,                       &
                  -one,                         &
                  wf%t1am(index1_first, 1),     & ! t_a_i
                  wf%n_v,                       &
                  x_ib_cd,                      &
                  length_2*length_4*length_3,   &
                  one,                          &
                  g_vv_vv,                      &
                  length_1)
!
      call deallocator(x_ib_cd, wf%n_o*length_2, length_3*length_4)
!
!     :: Term 2 and 3 ::
!     - sum_(J) sum_(k) t_c_k * L_ab_J * L_kd_J
!       sum_(ik) t_a_i * t_c_k * (sum_(J) L_ib_J * L_kd_J) 
!
      call allocator(L_dk_J, wf%n_o*length_4, wf%n_J)
      call allocator(L_ab_J, length_1*length_2, wf%n_J)
!
      call wf%read_cholesky_ai(L_dk_J, index4_first, index4_last, 1, wf%n_o)
      call wf%read_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
!
      call allocator(x_ab_dk, length_1*length_2, wf%n_o*length_4)
!
!     x_ab_dk = sum_(J)L_ab_J * L_dk_J
!
      call dgemm('N', 'T',             &
                  length_1*length_2,   &
                  wf%n_o*length_4,     &
                  wf%n_J,              &
                  one,                 &
                  L_ab_J,              &
                  length_1*length_2,   &
                  L_dk_J,              &
                  wf%n_o*length_4,     &
                  zero,                &
                  x_ab_dk,             &
                  length_1*length_2)
!
      call deallocator(L_dk_J, wf%n_o*length_2, wf%n_J)
!
!     g_ab_dc = sum_(k)t_c_k * x_ab_dk
!
      call allocator(g_ab_dc, length_1*length_2, length_4*length_3)
!
      call dgemm('N', 'T',                   &
                  length_1*length_2*length_4,&
                  length_3,                  &
                  wf%n_o,                    &
                  -one,                      &
                  x_ab_dk,                   &
                  length_1*length_2*length_4,&
                  wf%t1am(index3_first, 1),  & ! t_c_k
                  wf%n_v,                    &
                  zero,                      &
                  g_ab_dc,                   &
                  length_1*length_2*length_4)
!
      call deallocator(x_ab_dk, length_1*length_2, wf%n_o*length_4)
!
      call allocator(x_ib_dk, wf%n_o*length_2, wf%n_o*length_4)
!
!     x_ib_dk = sum_(J)L_ib_J * L_dk_J
!
      call dgemm('N', 'T',          &
                  wf%n_o*length_2,  &
                  wf%n_o*length_4,  &
                  wf%n_J,           &
                  one,              &
                  L_ib_J,           &
                  wf%n_o*length_2,  &
                  L_dk_J,           &
                  wf%n_o*length_4,  &
                  zero,             &
                  x_ib_dk,          &
                  wf%n_o*length_2)
!
      call deallocator(L_dk_J, wf%n_o*length_4, wf%n_J)
      call deallocator(L_ib_J, wf%n_o*length_2, wf%n_J)
!
      call allocator(x_ib_dc, wf%n_o*length_2, length_4*length_3)
!
!     x_ib_dc = sum_(k)t_c_k * x_ib_dk
!
      call dgemm('N', 'T',                   &
                  wf%n_o*length_2*length_4,  &
                  length_3,                  &
                  wf%n_o,                    &
                  one,                       &
                  x_ib_dk,                   &
                  wf%n_o*length_2*length_4,  &
                  wf%t1am(index3_first, 1),  & ! t_c_k
                  wf%n_v,                    &
                  zero,                      &
                  x_ib_dc,                   &
                  wf%n_o*length_2*length_4)
!
      call deallocator(x_ib_dk, wf%n_o*length_2, wf%n_o*length_4)
!
!     g_ab_dc += sum_(i) t_a_i * x_ib_dc
!
      call dgemm('N', 'N',                      &
                  length_1,                     &
                  length_2*length_4*length_3,   &
                  wf%n_o,                       &
                  one,                          &
                  wf%t1am(index1_first, 1),     & ! t_a_i
                  wf%n_v,                       &
                  x_ib_dc,                      &
                  length_2*length_4*length_3,   &
                  one,                          &
                  g_ab_dc,                      &
                  length_1)
!
      call deallocator(x_ib_dc, wf%n_o*length_2, length_4*length_3)
!
!     g_vv_vv = g_ab_cd(ab, cd) += g_ab_dc(ab, dc)
!
      do c = 1, length_3
         do d = 1, length_4
!
            cd = index_two(c, d, length_3)
            dc = index_two(d, c, length_4)
!
            g_vv_vv(:,cd) = g_vv_vv(:,cd) + g_ab_dc(:, dc)
!
         enddo
      enddo
!
      call deallocator(g_ab_dc, length_1*length_2, length_4*length_3)
!
   end subroutine t1_transform_vv_vv_ccs
!
   module subroutine store_electronic_repulsion_integrals_ccs(wf)
!!
!!    Store Electronic Repulsion Integrals 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Tests whether it is possible to store vir-vir-vir-vir integrals and,
!!    if possible, writes the integrals to disk 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      integer  :: required_space 
      real(dp) :: required_space_gb
!
      integer(i15) :: unit_g_abcd = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1     ! Error integer for file handling
      integer(i15) :: rec_number = -1  ! The record where g_abcd is positioned 
!
!     Batching variables
!
      integer(i15) :: required_mem = -1
      integer(i15) :: available_mem = -1 
!
      integer(i15) :: b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0, b_n_batch = 0
      integer(i15) :: d_first = 0, d_last = 0, d_length = 0, d_max_length = 0, d_batch = 0
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0, bcd = 0, cd_packed = 0, I = 0
!
      real(dp) :: begin_timer, end_timer
!
!     Cholesky vectors 
!
      real(dp), dimension(:,:), allocatable :: L_ab_J 
      real(dp), dimension(:,:), allocatable :: L_cd_J
!
!     The electronic repulsion integral 
!
      real(dp), dimension(:,:), allocatable :: g_a_bcd ! g_abcd 
!
!     Calculate the disk space (in GB) required to store the vir-vir-vir-vir integrals 
!
!     For non-T1-transformed integrals, there is a four-fold symmetri, 
!
!        ((a >= b) >= (c >= d),
!
!     giving required space = (n_v*(n_v+1)/2)*(n_v*(n_v+1)/2 + 1)/2
!
      required_space = ((wf%n_v**3)*(wf%n_v+1))/2 
!
!     This is the required space in number of double precision numbers (8 bytes per such number).
!     We convert this number to gigabytes. 
!
      required_space = 8*required_space ! in bytes
!
      required_space_gb = real(required_space)*(1.0D-9)
!
!     Test whether there is room for the integrals & save if this is the case 
!
      if (required_space_gb .lt. wf%settings%disk_space) then 
!
         call cpu_time(begin_timer)
!
!        Open file for writing integrals - one record: (a, bcd) = (1:n_v, bcd)
!
         call generate_unit_identifier(unit_g_abcd)
         open(unit=unit_g_abcd, file='g_abcd', action='write', status='unknown', &
               access='direct', form='unformatted', recl=dp*(wf%n_v), iostat=ioerror)
!
!        In calculating g_ab_cd, we will batch over the b and d indices 
!
         required_mem = max(2*(wf%n_v)**2*(wf%n_J) + 4*(wf%n_v)*(wf%n_o)*(wf%n_J), & ! Needed to get L_ab_J & L_cd_J
                           (wf%n_v)**4 + 2*(wf%n_v)**2*(wf%n_J)) ! Needed to get g_ac_bd
!         
         required_mem  = 4*required_mem
         available_mem = get_available()
!
         b_max_length = 0
         call num_two_batch(required_mem, available_mem, b_max_length, b_n_batch, wf%n_v)
!
         do b_batch = 1, b_n_batch ! Batch over b index 
!
            call batch_limits(b_first, b_last, b_batch, b_max_length, wf%n_v)
            b_length = b_last - b_first + 1  
!
!           Start looping over batches of d
!
            d_first  = 0
            d_last   = 0
            d_length = 0
!
            d_max_length = b_max_length
!
            call allocator(L_ab_J, (wf%n_v)*b_length, wf%n_J)
            call wf%read_cholesky_ab(L_ab_J, 1, wf%n_v, b_first, b_last)
!
            do d_batch = 1, b_batch ! Batch over d index; restricted by b index 
!
               call batch_limits(d_first, d_last, d_batch, d_max_length, wf%n_v)
               d_length = d_last - d_first + 1 
!
!              Calculate the integrals g_ab_cd, where b and d are restricted by
!              the current batching limits 
!
               call allocator(L_cd_J, (wf%n_v)*d_length, wf%n_J)
!
               call wf%read_cholesky_ab(L_cd_J, 1, wf%n_v, d_first, d_last)
!
               call allocator(g_a_bcd, (wf%n_v), b_length*(wf%n_v)*d_length)
!
               call dgemm('N', 'T',           &
                           (wf%n_v)*b_length, & 
                           (wf%n_v)*d_length, &
                           wf%n_J,            &
                           one,               &
                           L_ab_J,            &
                           (wf%n_v)*b_length, &
                           L_cd_J,            &
                           (wf%n_v)*d_length, &
                           zero,              &
                           g_a_bcd,           & ! g_ab_cd 
                           (wf%n_v)*b_length)
!
               call deallocator(L_cd_J, (wf%n_v)*d_length, wf%n_J)
!
!              Save the integrals to disk 
!
               do d = 1, d_length
                  do c = 1, wf%n_v
                     do b = 1, b_length
!
!                       Calculate record number 
!
                        cd_packed = index_packed(c, d + d_first -1)
                        bcd = index_two(b + b_first - 1, cd_packed, wf%n_v)
!
                        if (c .ge. (d + d_first - 1)) then 
!
!                          Write integrals to that record 
!
                           write(unit_g_abcd, rec=bcd) (g_a_bcd(I, bcd), I = 1, wf%n_v)
!
                        endif
!
                     enddo
                  enddo
               enddo
!
               call deallocator(g_a_bcd, (wf%n_v), b_length*(wf%n_v)*d_length)
!
            enddo ! End of batches over d
!
            call deallocator(L_ab_J, (wf%n_v)*b_length, wf%n_J)
!
         enddo ! End of batches over b 
!
!        Test for file handling error 
!
         if (ioerror .ne. 0) write(unit_output,'(t3,a)') 'Error: write error in store_electronic_repulsion_integrals_ccs'
!
!        Close file 
!
         close(unit_g_abcd)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then
! 
            write(unit_output,'(t3,a36,f14.8)') 'Time used to store g_abcd (seconds):', end_timer - begin_timer
            flush(unit_output)
!
         endif
!
      endif 
!
   end subroutine store_electronic_repulsion_integrals_ccs
!
!
end submodule integrals