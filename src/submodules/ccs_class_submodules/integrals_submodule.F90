submodule (ccs_class) integrals
!
!!
!!    Integrals submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
!!    Contains procedures for construction of general integrals, and 
!!    specifically:
!!
!!    - electronic repulsion integrals
!!    - (...no other integrals yet...)
!!    
!!    o - occupied index
!!    v - virtual index
!!
!!    Note: for normal use, the get_pq_rs routines should be called. These
!!    will call the appropriate routines for constructing (or reading) the integrals.
!!    For example, 
!!
!!       integral_type = 'electronic_repulsion'
!!       call wf%get_vo_vo(integral_type, g_ai_bj)
!!
!!    will place the T1-transformed g_aibj integrals in the g_ai_bj array.
!!    The indices may also be restricted (though this is optional):
!!
!!       integral_type = 'electronic_repulsion'
!!       call wf%get_vo_vo(integral_type, g_ai_bj, a_first, a_last, i_first, i_last, ...)
!!
!
contains
!
!    -::- Get integral routines -::-
!    :::::::::::::::::::::::::::::::
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
      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
               .and. present(index2_first) .and. present(index2_last) &
               .and. present(index3_first) .and. present(index3_last) &
               .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_o
         local_index3_last = wf%n_o
         local_index4_last = wf%n_o
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_oo_oo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_oo_electronic_repulsion(x_oo_oo,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_oo_oo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_o
         local_index3_last = wf%n_o
         local_index4_last = wf%n_v
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_oo_ov'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_ov_electronic_repulsion(x_oo_ov,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_oo_ov'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_v
         local_index3_last = wf%n_o
         local_index4_last = wf%n_o
!
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_ov_oo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_oo_electronic_repulsion(x_ov_oo,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_ov_oo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
               .and. present(index2_first) .and. present(index2_last) &
               .and. present(index3_first) .and. present(index3_last) &
               .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_o
         local_index3_last = wf%n_v
         local_index4_last = wf%n_o
!         
      else
         write(unit_output,*) 'Error: some optionals missing in call to get_oo_vo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_vo_electronic_repulsion(x_oo_vo,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_oo_vo'
         stop
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
      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_o
         local_index3_last = wf%n_o
         local_index4_last = wf%n_o
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_vo_oo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
        call wf%get_vo_oo_electronic_repulsion(x_vo_oo,                      & 
                                      local_index1_first, local_index1_last, &
                                      local_index2_first, local_index2_last, &
                                      local_index3_first, local_index3_last, &
                                      local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vo_oo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_o
         local_index3_last = wf%n_v
         local_index4_last = wf%n_v
!         
      else
         write(unit_output,*) 'Error: some optionals missing in call to get_oo_vv'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_oo_vv_electronic_repulsion(x_oo_vv,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_oo_vv'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_v
         local_index3_last = wf%n_o
         local_index4_last = wf%n_o
!         
      else
         write(unit_output,*)'WARNING: Some optionals missing in call to get_vv_oo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_oo_electronic_repulsion(x_vv_oo,          & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'WARNING: unknown integral type requested from get_vv_oo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_v
         local_index3_last = wf%n_o
         local_index4_last = wf%n_v
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_ov_ov'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_ov_electronic_repulsion(x_ov_ov,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_ov_ov'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_o
         local_index3_last = wf%n_v
         local_index4_last = wf%n_o
!         
      else
!
         write(unit_output,*) 'Error: Some optionals missing in call to get_vo_vo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vo_vo_electronic_repulsion(x_vo_vo,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vo_vo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_v
         local_index3_last = wf%n_v
         local_index4_last = wf%n_o
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_ov_vo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_vo_electronic_repulsion(x_ov_vo,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_ov_vo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
               .and. present(index2_first) .and. present(index2_last) &
               .and. present(index3_first) .and. present(index3_last) &
               .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_o
         local_index3_last = wf%n_o
         local_index4_last = wf%n_v
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_vo_ov'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vo_ov_electronic_repulsion(x_vo_ov,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vo_ov'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_o
         local_index2_last = wf%n_v
         local_index3_last = wf%n_v
         local_index4_last = wf%n_v
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_ov_vv'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_ov_vv_electronic_repulsion(x_ov_vv,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_ov_vv'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_v
         local_index3_last = wf%n_o
         local_index4_last = wf%n_v
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_vv_ov'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_ov_electronic_repulsion(x_vv_ov,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vv_ov'
         stop
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
      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_o
         local_index3_last = wf%n_v
         local_index4_last = wf%n_v
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_vo_vv'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vo_vv_electronic_repulsion(x_vo_vv,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vo_vv'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_v
         local_index3_last = wf%n_v
         local_index4_last = wf%n_o
!         
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_vv_vo'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_vo_electronic_repulsion(x_vv_vo,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vv_vo'
         stop
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

      integer(i15) :: local_index1_first, local_index1_last
      integer(i15) :: local_index2_first, local_index2_last
      integer(i15) :: local_index3_first, local_index3_last
      integer(i15) :: local_index4_first, local_index4_last
!
!     Set local index variables.
!     Necessary because optional arguments cannot have default values
!
      if (     present(index1_first) .and. present(index1_last) &
         .and. present(index2_first) .and. present(index2_last) &
         .and. present(index3_first) .and. present(index3_last) &
         .and. present(index4_first) .and. present(index4_last) ) then
!
         local_index1_first = index1_first
         local_index2_first = index2_first
         local_index3_first = index3_first
         local_index4_first = index4_first
!
         local_index1_last = index1_last
         local_index2_last = index2_last
         local_index3_last = index3_last
         local_index4_last = index4_last
!

      elseif ( .not. (present(index1_first) .and. present(index1_last) &
                .and. present(index2_first) .and. present(index2_last) &
                .and. present(index3_first) .and. present(index3_last) &
                .and. present(index4_first) .and. present(index4_last) )) then
!
         local_index1_first = 1
         local_index2_first = 1
         local_index3_first = 1
         local_index4_first = 1
!
         local_index1_last = wf%n_v
         local_index2_last = wf%n_v
         local_index3_last = wf%n_v
         local_index4_last = wf%n_v
!
      else
!
         write(unit_output,*) 'Error: some optionals missing in call to get_vv_vv'
         stop
!
      endif
!
      if (trim(integral_type) == 'electronic_repulsion') then
!
         call wf%get_vv_vv_electronic_repulsion(x_vv_vv,                      & 
                                       local_index1_first, local_index1_last, &
                                       local_index2_first, local_index2_last, &
                                       local_index3_first, local_index3_last, &
                                       local_index4_first, local_index4_last)
!
      else
!
         write(unit_output,*) 'Error: unknown integral type requested from get_vv_vv'
         stop
!
      endif
!
   end subroutine get_vv_vv_ccs
!
!
!    -::- Get electronic repulsion integral routines -::-
!    ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     These routines should, normally, not be called directly. 
!     If "integral_type" equals "electronic_repulsion", the get_pq_rs
!     routines will call the routines below to get the correct integrals. 
!
!     Depending on type of calculation, integrals are either read from file 
!     or constructed directly from the Cholesky vectors.
!
!
   module subroutine get_vo_vo_electronic_repulsion_ccs(wf, x_vo_vo, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vo,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,vo (ordered as g_vo_vo)
!!
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_bj_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ai_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_bj_J, length_3*length_4, wf%n_J)
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
      call wf%mem%dealloc(L_ai_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_bj_J, length_3*length_4, wf%n_J)
!
   end subroutine get_vo_vo_electronic_repulsion_ccs
!
!
   module subroutine get_ov_vo_electronic_repulsion_ccs(wf, x_ov_vo, & 
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_bj_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_bj_J, length_3*length_4, wf%n_J)
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
      call wf%mem%dealloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_bj_J, length_3*length_4, wf%n_J)
!
   end subroutine get_ov_vo_electronic_repulsion_ccs
!
!
   module subroutine get_vo_ov_electronic_repulsion_ccs(wf, x_vo_ov, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vo,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,ov (ordered as g_vo_ov)
!!
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_jb_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      logical :: voov_t1_on_file
!
!     Test if we have T1-transformed on file
!
      inquire(file='g_t1_aijb',exist=voov_t1_on_file)
!
      if (voov_t1_on_file .and. wf%tasks%current .ne. 'ground_state') then
!  
         call wf%read_t1_vo_ov_electronic_repulsion(x_vo_ov,      & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
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
         call wf%mem%alloc(L_ai_J, length_1*length_2, wf%n_J)
         call wf%mem%alloc(L_jb_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ia(L_jb_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
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
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_ai_J, length_1*length_2, wf%n_J)
         call wf%mem%dealloc(L_jb_J, length_3*length_4, wf%n_J)
!
      endif
!
   end subroutine get_vo_ov_electronic_repulsion_ccs
!
!
   module subroutine get_ov_vv_electronic_repulsion_ccs(wf, x_ov_vv, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,vv (ordered as g_vo_vv)
!!
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_bc_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_bc_J, length_3*length_4, wf%n_J)
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
      call wf%mem%dealloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_bc_J, length_3*length_4, wf%n_J)
!
   end subroutine get_ov_vv_electronic_repulsion_ccs
!
!
   module subroutine get_vv_ov_electronic_repulsion_ccs(wf, x_vv_ov, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vv,ov integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,ov (ordered as g_vv_ov)
!!
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_ic_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      logical :: vvov_t1_on_file = .false.
!
!     Test if we have T1-transformed on file
!
      inquire(file='g_t1_bcia',exist=vvov_t1_on_file)
!
      if (vvov_t1_on_file .and. wf%tasks%current .ne. 'ground_state') then
!  
         call wf%read_t1_vv_ov_electronic_repulsion(x_vv_ov,      & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
      else
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
         call wf%mem%alloc(L_ab_J, length_1*length_2, wf%n_J)
         call wf%mem%alloc(L_ic_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ia(L_ic_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
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
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_ab_J, length_1*length_2, wf%n_J)
         call wf%mem%dealloc(L_ic_J, length_3*length_4, wf%n_J)
!
      endif
!
   end subroutine get_vv_ov_electronic_repulsion_ccs
!
!
   module subroutine get_vo_vv_electronic_repulsion_ccs(wf, x_vo_vv, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vo,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vo,vv (ordered as g_vo_vv)
!!
!!    Parameters
!! 
!!    index1_first, index1_last,   
!!    index2_first, index2_last,   
!!    index3_first, index3_last  and 
!!    index4_first, index4_last
!!
!!    are used to restrict indices of the integral.
!!
!!    Note: in an excited state calculation, the T1-transformed integrals 
!!    will be read from a file containing the vv_vo integrals. This involves
!!    both reading & reordering, thus requiring twice the memory: 2 * v^3 * o.
!!
!!    To minimize the number of wasteful reads from file, batching should be in 
!!    the third index, i.e. the capitalized letter in: g_vo_Vv.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: x_vo_vv
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_bc_J
!
      real(dp), dimension(:,:), allocatable :: x_vv_vo ! g_ab_ci, for reading from file
      integer(i15) :: I = 0, J = 0
!
      logical :: t1_vvvo_on_file = .false.
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Check if the t1-transformed integrals are on file 
!
      inquire(file='g_t1_abci',exist=t1_vvvo_on_file)
!
      if (t1_vvvo_on_file .and. wf%tasks%current .ne. 'ground_state') then 
!
!        Read the integrals from file 
!
         call wf%mem%alloc(x_vv_vo, length_3*length_4, length_1*length_2)
!
         call wf%read_t1_vv_vo_electronic_repulsion(x_vv_vo, &
                                 index3_first, index3_last,  &
                                 index4_first, index4_last,  &
                                 index1_first, index1_last,  &
                                 index2_first, index2_last)
!
!        Place the integrals in the correct positions 
!
         do I = 1, length_1*length_2
            do J = 1, length_3*length_4
!
               x_vo_vv(I,J) = x_vv_vo(J,I)
!
            enddo
         enddo
!
         call wf%mem%dealloc(x_vv_vo, length_3*length_4, length_1*length_2)
!
      else
!
!        Alllocate Cholesky vectors
!
         call wf%mem%alloc(L_ai_J, length_1*length_2, wf%n_J)
         call wf%mem%alloc(L_bc_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ab(L_bc_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
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
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_ai_J, length_1*length_2, wf%n_J)
         call wf%mem%dealloc(L_bc_J, length_3*length_4, wf%n_J)
!
      endif
!
   end subroutine get_vo_vv_electronic_repulsion_ccs
!
!
   module subroutine get_vv_vo_electronic_repulsion_ccs(wf, x_vv_vo, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vv,vo integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,vo (ordered as g_vv_vo)
!!
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_ci_J
!
      logical :: t1_vvvo_on_file = .false.
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
!     Lengths
!
      length_1 = index1_last - index1_first + 1
      length_2 = index2_last - index2_first + 1
      length_3 = index3_last - index3_first + 1
      length_4 = index4_last - index4_first + 1
!
!     Check if the t1-transformed integrals are on file 
!
      inquire(file='g_t1_abci',exist=t1_vvvo_on_file)
!
      if (t1_vvvo_on_file .and. wf%tasks%current .ne. 'ground_state') then 
!
!        Read the integrals from file 
!
         call wf%read_t1_vv_vo_electronic_repulsion(x_vv_vo, &
                                 index1_first, index1_last,  &
                                 index2_first, index2_last,  &
                                 index3_first, index3_last,  &
                                 index4_first, index4_last)
!
      else
!
!        Alllocate Cholesky vectors
!
         call wf%mem%alloc(L_ab_J, length_1*length_2, wf%n_J)
         call wf%mem%alloc(L_ci_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
         call wf%get_cholesky_ai(L_ci_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
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
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_ab_J, length_1*length_2, wf%n_J)
         call wf%mem%dealloc(L_ci_J, length_3*length_4, wf%n_J)
!
      endif
!
   end subroutine get_vv_vo_electronic_repulsion_ccs
!
!
!
   module subroutine get_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv, & 
                                       index1_first, index1_last,    &
                                       index2_first, index2_last,    &
                                       index3_first, index3_last,    &
                                       index4_first, index4_last)
!!
!!    Get g_vv,vv integral (CCS),
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep. 2017.
!! 
!!    Construct two-electron repulsion integral g_vv,vv (ordered as g_vv_vv)
!!
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_cd_J
!
      integer(i15) :: ab = 0, cd = 0
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      logical :: vvvv_on_file    = .false.
      logical :: t1_vvvv_on_file = .false.
!
      real(dp) :: begin_timer, end_timer
!
      inquire(file='g_abcd',exist=vvvv_on_file)
      inquire(file='g_t1_abcd',exist=t1_vvvv_on_file)
!
      if (vvvv_on_file) then
!
!        Read x_vv_vv
!
         call cpu_time(begin_timer)
!
         call wf%read_vv_vv_electronic_repulsion(x_vv_vv, &
                                 index1_first, index1_last, &
                                 index2_first, index2_last, &
                                 index3_first, index3_last, &
                                 index4_first, index4_last)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then 
!
   !         write(unit_output,'(t6,a27,f14.8)') 'Read abcd (seconds):', end_timer-begin_timer
!
         endif
!
!        T1-transform x_vv_vv
!
         call cpu_time(begin_timer)
!
         call wf%t1_transform_vv_vv(x_vv_vv,                &
                                 index1_first, index1_last, &
                                 index2_first, index2_last, &
                                 index3_first, index3_last, &
                                 index4_first, index4_last)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then 
!
   !         write(unit_output,'(t6,a27,f14.8)') 't1-transform (seconds):', end_timer-begin_timer
!
         endif 
!
      elseif (t1_vvvv_on_file .and. wf%tasks%current .ne. 'ground_state') then 
!
!        Read x_vv_vv
!
         call cpu_time(begin_timer)
!
         call wf%read_t1_vv_vv_electronic_repulsion(x_vv_vv, &
                                 index1_first, index1_last,  &
                                 index2_first, index2_last,  &
                                 index3_first, index3_last,  &
                                 index4_first, index4_last)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then 
!
   !         write(unit_output,'(t6,a27,f14.8)') 'Read t1 abcd (seconds):', end_timer-begin_timer 
!
         endif        
!
      else
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
         call wf%mem%alloc(L_ab_J, length_1*length_2, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
!
!        Alllocate Cholesky vectors
!
         call wf%mem%alloc(L_cd_J, length_3*length_4, wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_cd_J, index3_first, index3_last, index4_first, index4_last)
!
!        Construct integral
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
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_ab_J, length_1*length_2, wf%n_J)
         call wf%mem%dealloc(L_cd_J, length_3*length_4, wf%n_J)
!
      endif
!
   end subroutine get_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine read_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv, & 
                                       index1_first, index1_last,     &
                                       index2_first, index2_last,     &
                                       index3_first, index3_last,     &
                                       index4_first, index4_last)
!!
!!    Read vvvv Electronic Repulsion (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Reads the non-T1-transformed g_abcd integrals from file,
!!    with indices a,b,c and d restricted as requested. 
!!
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
!     Index lengths 
!
      integer(i15) :: length_a = 0, length_b = 0, length_c = 0, length_d = 0
!
!     Calculate lengths of indices a,b,c,d 
!
      length_a = index1_last - index1_first + 1
      length_b = index2_last - index2_first + 1
      length_c = index3_last - index3_first + 1
      length_d = index4_last - index4_first + 1
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
      call wf%mem%alloc(x_v, wf%n_v, 1)
      x_v = zero 
!
      if (length_a .ne. wf%n_v) then ! Batching over first index, a 
!
         if (length_b .ne. wf%n_v) then ! Batching over second index, b, as well
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
                     x_v = zero
                     read(unit_g_abcd, rec=bcd, iostat=ioerror) (x_v(a, 1), a = 1, wf%n_v)
!
!                    Place the integral into x_vv_vv = g_ab_cd 
!
                     cd = index_two(c - index3_first + 1, &
                                    d - index4_first + 1, &
                                    index3_last - index3_first + 1)
!
                     do a = index1_first, index1_last
!
                        ab = index_two(a - index1_first + 1, &
                                       b - index2_first + 1, &
                                       index1_last - index1_first + 1)
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
                     read(unit_g_abcd, rec=acd, iostat=ioerror) (x_v(b, 1), b = 1, wf%n_v)
!
!                    Place the integral into x_vv_vv = g_ab_cd 
!
                     cd = index_two(c - index3_first + 1, &
                                    d - index4_first + 1, &
                                    index3_last - index3_first + 1)
!
                     do b = index2_first, index2_last
!
                        ab = index_two(a - index1_first + 1, &
                                       b - index2_first + 1, &
                                       index1_last - index1_first + 1)
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
      else ! No batching over first index
!
         do d = 1, length_d
            do c = 1, length_c
               do b = 1, length_b
!
                  bcd = index_two(b + index2_first - 1, index_packed(c,d), wf%n_v) ! Record number 
!
!                 Read g_a_bcd, a = 1, n_v, into x_v 
!
                  read(unit_g_abcd, rec=bcd, iostat=ioerror) (x_v(a, 1), a = 1, wf%n_v)
!
!                 Place the integral into x_vv_vv = g_ab_cd 
!
                  cd = index_two(c, d, length_c)
!
                  do a = 1, length_a
!
                     ab = index_two(a, b, length_a)
!
                     x_vv_vv(ab, cd) = x_v(a + index1_first - 1, 1)
!
                  enddo
!
               enddo
            enddo
         enddo 
!
      endif 
!
!     Deallocate temporary vector 
!
      call wf%mem%dealloc(x_v, wf%n_v, 1)
!
!     Close file containing the g_abcd integrals 
!
      close(unit_g_abcd)
!
   end subroutine read_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine read_t1_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv, & 
                                       index1_first, index1_last,        &
                                       index2_first, index2_last,        &
                                       index3_first, index3_last,        &
                                       index4_first, index4_last)
!!
!!    Read T1 vvvv Electronic Repulsion (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Reads the T1-transformed g_abcd integrals from file,
!!    with indices a,b,c and d restricted as requested. 
!!
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
      integer(i15) :: unit_g_t1_abcd = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1     ! Error integer for file handling
      integer(i15) :: rec_number = -1  ! The record where g_abcd is positioned 
!
!     Indices 
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0, ab = 0, cd = 0, bcd = 0, acd = 0 
!
!     Index lengths 
!
      integer(i15) :: length_a = 0, length_b = 0, length_c = 0, length_d = 0
!
!     Calculate lengths of indices a,b,c,d 
!
      length_a = index1_last - index1_first + 1
      length_b = index2_last - index2_first + 1
      length_c = index3_last - index3_first + 1
      length_d = index4_last - index4_first + 1
!
!     Open file containing the g_abcd integrals, ordered as
!     g_a_bcd, where bcd index is unpakced 
!
!     The compound index bcd determines the record number,
!     where the record includes the integrals g(a,bcd), a = 1, n_v
!
      call generate_unit_identifier(unit_g_t1_abcd)
      open(unit=unit_g_t1_abcd, file='g_t1_abcd', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_v), iostat=ioerror)
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,*) 'Error: could not open file g_t1_abcd in read_t1_vv_vv_electronic_repulsion_ccs'
         stop
!
      endif
!
!     Allocate the integral x_v the integrals g_a_bcd for a given bcd, a = 1, n_v 
!
      call wf%mem%alloc(x_v, wf%n_v, 1)
      x_v = zero 
!
      if (length_a .ne. wf%n_v) then ! Batching over first index, a 
!
         if (length_b .ne. wf%n_v) then ! Batching over second index, b, as well
!
!           No simple tricks available. Read all integrals g_a_bcd, a = 1, n_v, into x_v, then through away!
!
            do d = index4_first, index4_last
               do c = index3_first, index3_last
                  do b = index2_first, index2_last
!
                     bcd = index_three(b, c, d, wf%n_v, wf%n_v) ! Record number 
!
!                    Read g_a_bcd, a = 1, n_v, into x_v 
!
                     x_v = zero
                     read(unit_g_t1_abcd, rec=bcd, iostat=ioerror) (x_v(a, 1), a = 1, wf%n_v)
!
!                    Place the integral into x_vv_vv = g_ab_cd 
!
                     cd = index_two(c - index3_first + 1, &
                                    d - index4_first + 1, &
                                    index3_last - index3_first + 1)
!
                     do a = index1_first, index1_last
!
                        ab = index_two(a - index1_first + 1, &
                                       b - index2_first + 1, &
                                       index1_last - index1_first + 1)
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
                     acd = index_three(a, c, d, wf%n_v, wf%n_v) ! Record number 
!
!                    Read g_b_acd, b = 1, n_v, into x_v 
!
                     read(unit_g_t1_abcd, rec=acd, iostat=ioerror) (x_v(b, 1), b = 1, wf%n_v)
!
!                    Place the integral into x_vv_vv = g_ab_cd 
!
                     cd = index_two(c - index3_first + 1, &
                                    d - index4_first + 1, &
                                    index3_last - index3_first + 1)
!
                     do b = index2_first, index2_last
!
                        ab = index_two(a - index1_first + 1, &
                                       b - index2_first + 1, &
                                       index1_last - index1_first + 1)
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
      else ! No batching over first index
!
         do d = 1, length_d
            do c = 1, length_c
               do b = 1, length_b
!
                  bcd = index_three(b + index2_first - 1, c, d, wf%n_v, wf%n_v) ! Record number 
!
!                 Read g_a_bcd, a = 1, n_v, into x_v 
!
                  read(unit_g_t1_abcd, rec=bcd, iostat=ioerror) (x_v(a, 1), a = 1, wf%n_v)
!
!                 Place the integral into x_vv_vv = g_ab_cd 
!
                  cd = index_two(c, d, length_c)
!
                  do a = 1, length_a
!
                     ab = index_two(a, b, length_a)
!
                     x_vv_vv(ab, cd) = x_v(a + index1_first - 1, 1)
!
                  enddo
!
               enddo
            enddo
         enddo 
!
      endif 
!
!     Deallocate temporary vector 
!
      call wf%mem%dealloc(x_v, wf%n_v, 1)
!
!     Close file containing the g_abcd integrals 
!
      close(unit_g_t1_abcd)
!
   end subroutine read_t1_vv_vv_electronic_repulsion_ccs
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_kl_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_kl_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ij(L_kl_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_kl_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_ka_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
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
      call wf%mem%alloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_ka_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ia(L_ka_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_ka_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_jk_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_jk_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ia(L_ia_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ij(L_jk_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_jk_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_ak_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_ak_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ai(L_ak_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_ak_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ai_J, L_jk_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ai_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_jk_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ai(L_ai_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ij(L_jk_J)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ai_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_jk_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ij_J, L_ab_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_ab_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ij(L_ij_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ab(L_ab_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ij_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_ab_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ab_J, L_ij_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ab_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_ij_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ij(L_ij_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ab_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_ij_J, length_3*length_4, wf%n_J)
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
!!    Parameters
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
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp), dimension(:,:), allocatable :: L_ia_J, L_jb_J
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
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
      call wf%mem%alloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%alloc(L_jb_J, length_3*length_4, wf%n_J)
!
!     Get T1-transformed Cholesky vectors
!
      call wf%get_cholesky_ia(L_ia_J, index1_first, index1_last, index2_first, index2_last)
      call wf%get_cholesky_ia(L_jb_J, index3_first, index3_last, index4_first, index4_last)
!
!     Construct integral
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
!     Deallocate Cholesky vectors
!
      call wf%mem%dealloc(L_ia_J, length_1*length_2, wf%n_J)
      call wf%mem%dealloc(L_jb_J, length_3*length_4, wf%n_J)
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
!!       g_ab_cd_T1 = g_ab_cd - sum_(J) sum_(i) t_a_i * L_ib_J * L_cd_J
!!                            - sum_(J) sum_(k) t_c_k * L_kd_J * L_ab_J
!!                            + sum_(J) sum_(ki) t_a_i * t_c_k * L_kd_J * L_ib_J
!! 
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: g_vv_vv
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      integer(i15) :: a = 0, b = 0,c = 0, d = 0, ab = 0, cd = 0, dc = 0, i = 0, j = 0
!
      integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      real(dp), dimension(:,:), allocatable :: L_ib_J, L_dk_J, L_cd_J, L_ab_J
      real(dp), dimension(:,:), allocatable :: x_ib_dk, x_ib_dc, x_ib_cd, x_ab_dk
      real(dp), dimension(:,:), allocatable :: g_ab_dc, g_ab_cd
!
      length_1 = index1_last - index1_first + 1 ! a 
      length_2 = index2_last - index2_first + 1 ! b 
      length_3 = index3_last - index3_first + 1 ! c
      length_4 = index4_last - index4_first + 1 ! d
!
!        Test if we are batching
!
         if ((length_1 .eq. wf%n_v) .and. (length_2 .eq. wf%n_v) .and. &
            (length_3 .eq. wf%n_v) .and. (length_4 .eq. wf%n_v) ) then
!
!        We are NOT batching
!
!        :: Term 1 and 2 ::
!
!         - sum_(J) sum_(i) t_a_i * L_ib_J * L_cd_J - sum_(J) sum_(k) t_c_k * L_ab_J * L_kd_J
!        = - sum_(J) sum_(i) (t_a_i * L_ib_J * L_cd_J -  t_c_i * L_ab_J * L_id_J)
!        = - sum_(i) (t_a_i * x_ib_cd -  t_c_i * x_id_ab)
!        = - g_ab_cd - g_cd_ab
!
!
         call wf%mem%alloc(L_ib_J, (wf%n_o)*length_2, wf%n_J)
         call wf%mem%alloc(L_cd_J, length_3*length_4, wf%n_J)
!
         call wf%read_cholesky_ia(L_ib_J, 1, wf%n_o, index2_first, index2_last)
         call wf%read_cholesky_ab(L_cd_J, index3_first, index3_last, index4_first, index4_last)
!
         call wf%mem%alloc(x_ib_cd, wf%n_o*length_2, length_3*length_4)
!
!        x_ib_cd = sum_(J) L_ib_J * L_cd_J
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
         call wf%mem%dealloc(L_cd_J, length_3*length_4, wf%n_J)
!
!        g_ab_cd -= sum_(i)t_a_i* x_ib_cd
!
         call wf%mem%alloc(g_ab_cd, length_1*length_2, length_3*length_4)
         call dgemm('N', 'N',                      &
                     length_1,                     &
                     length_2*length_4*length_3,   &
                     wf%n_o,                       &
                     -one,                         &
                     wf%t1am(index1_first, 1),     & ! t_a_i
                     wf%n_v,                       &
                     x_ib_cd,                      &
                     wf%n_o,                       &
                     zero,                         &
                     g_ab_cd,                      & ! g_a_bcd
                     length_1)
!
         call wf%mem%dealloc(x_ib_cd, wf%n_o*length_2, length_3*length_4)
!
!        Add to g_vv_vv
!        g_vv_vv = - g_ab_cd(ab, cd) - g_ab_cd(cd, ab)
!
!$omp parallel do schedule(static) private(a,c,d,ab,cd)
         do b = 1, length_2
            do a = 1, length_1
               ab = index_two(a, b, length_1)
               do d = 1, length_4
                  do c = 1, length_3
                     cd = index_two(c, d, length_3)
                     g_vv_vv(ab, cd) = g_vv_vv(ab, cd) + g_ab_cd(ab, cd) + g_ab_cd(cd, ab)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call wf%mem%dealloc(g_ab_cd, length_1*length_2, length_3*length_4)
!
!        :: Term 3 ::
!
!          sum_(ik) t_a_i * t_c_k * (sum_(J) L_ib_J * L_kd_J) 
!
         call wf%mem%alloc(L_dk_J, wf%n_o*length_4, wf%n_J)
!
         call wf%read_cholesky_ai(L_dk_J, index4_first, index4_last, 1, wf%n_o)
         call wf%mem%alloc(x_ib_dk, wf%n_o*length_2, wf%n_o*length_4)
!
!        x_ib_dk = sum_(J)L_ib_J * L_dk_J
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
         call wf%mem%dealloc(L_dk_J, wf%n_o*length_4, wf%n_J)
         call wf%mem%dealloc(L_ib_J, wf%n_o*length_2, wf%n_J)
!
         call wf%mem%alloc(x_ib_dc, wf%n_o*length_2, length_4*length_3)
!
!        x_ib_dc = sum_(k)t_c_k * x_ib_dk
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
         call wf%mem%dealloc(x_ib_dk, wf%n_o*length_2, wf%n_o*length_4)
!
!        g_ab_dc += sum_(i) t_a_i * x_ib_dc
!
         call wf%mem%alloc(g_ab_dc, length_1*length_2, length_4*length_3)
!
         call dgemm('N', 'N',                      &
                     length_1,                     &
                     length_2*length_4*length_3,   &
                     wf%n_o,                       &
                     one,                          &
                     wf%t1am(index1_first, 1),     & ! t_a_i
                     wf%n_v,                       &
                     x_ib_dc,                      &
                     wf%n_o,                       &
                     zero,                         & ! g_a_bdc
                     g_ab_dc,                      &
                     length_1)
!
         call wf%mem%dealloc(x_ib_dc, wf%n_o*length_2, length_4*length_3)
!
!        g_vv_vv = g_ab_cd(ab, cd) += g_ab_dc(ab, dc)
!
!$omp parallel do schedule(static) private(d,cd,dc)
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
!$omp end parallel do 
!
         call wf%mem%dealloc(g_ab_dc, length_1*length_2, length_4*length_3)
!
      else
!
!        We are batching
!
!        ::Term 1:: 
!        - sum_(J) sum_(i) t_a_i * L_ib_J * L_cd_J 
!
         call wf%mem%alloc(L_ib_J, (wf%n_o)*length_2, wf%n_J)
         call wf%mem%alloc(L_cd_J, length_3*length_4, wf%n_J)
!
         call wf%read_cholesky_ia(L_ib_J, 1, wf%n_o, index2_first, index2_last)
         call wf%read_cholesky_ab(L_cd_J, index3_first, index3_last, index4_first, index4_last)
!
         call wf%mem%alloc(x_ib_cd, wf%n_o*length_2, length_3*length_4)
!
!        x_ib_cd = sum_(J) L_ib_J * L_cd_J
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
         call wf%mem%dealloc(L_cd_J, length_3*length_4, wf%n_J)
         call wf%mem%dealloc(L_ib_J, wf%n_o*length_2, wf%n_J)
!
!        g_vv_vv = g_ab_cd -= sum_(i)t_a_i* x_ib_cd
!
         call dgemm('N', 'N',                      &
                     length_1,                     &
                     length_2*length_4*length_3,   &
                     wf%n_o,                       &
                     -one,                         &
                     wf%t1am(index1_first, 1),     & ! t_a_i
                     wf%n_v,                       &
                     x_ib_cd,                      &
                     wf%n_o,                       &
                     one,                          &
                     g_vv_vv,                      & ! g_a_bcd
                     length_1)
!
         call wf%mem%dealloc(x_ib_cd, wf%n_o*length_2, length_3*length_4)
!
!        :: Term 2 and 3 ::
!        - sum_(J) sum_(k) t_c_k * L_ab_J * L_kd_J
!          sum_(ik) t_a_i * t_c_k * (sum_(J) L_ib_J * L_kd_J) 
!
         call wf%mem%alloc(L_dk_J, wf%n_o*length_4, wf%n_J)
         call wf%mem%alloc(L_ab_J, length_1*length_2, wf%n_J)
!
         call wf%read_cholesky_ai(L_dk_J, index4_first, index4_last, 1, wf%n_o)
         call wf%read_cholesky_ab(L_ab_J, index1_first, index1_last, index2_first, index2_last)
!
         call wf%mem%alloc(x_ab_dk, length_1*length_2, wf%n_o*length_4)
!
!        x_ab_dk = sum_(J)L_ab_J * L_dk_J
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
         call wf%mem%dealloc(L_ab_J, length_1*length_2, wf%n_J)
!
!        g_ab_dc = - sum_(k)t_c_k * x_ab_dk
!
         call wf%mem%alloc(g_ab_dc, length_1*length_2, length_4*length_3)
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
         call wf%mem%dealloc(x_ab_dk, length_1*length_2, wf%n_o*length_4)
!
         call wf%mem%alloc(x_ib_dk, wf%n_o*length_2, wf%n_o*length_4)
!
!        x_ib_dk = sum_(J)L_ib_J * L_dk_J
!    
         call wf%mem%alloc(L_ib_J, (wf%n_o)*length_2, wf%n_J)
         call wf%read_cholesky_ia(L_ib_J, 1, wf%n_o, index2_first, index2_last)
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
         call wf%mem%dealloc(L_dk_J, wf%n_o*length_4, wf%n_J)
         call wf%mem%dealloc(L_ib_J, wf%n_o*length_2, wf%n_J)
!
         call wf%mem%alloc(x_ib_dc, wf%n_o*length_2, length_4*length_3)
!
!        x_ib_dc = sum_(k)t_c_k * x_ib_dk
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
         call wf%mem%dealloc(x_ib_dk, wf%n_o*length_2, wf%n_o*length_4)
!
!        g_ab_dc += sum_(i) t_a_i * x_ib_dc
!
         call dgemm('N', 'N',                      &
                     length_1,                     &
                     length_2*length_4*length_3,   &
                     wf%n_o,                       &
                     one,                          &
                     wf%t1am(index1_first, 1),     & ! t_a_i
                     wf%n_v,                       &
                     x_ib_dc,                      &
                     wf%n_o,                       &
                     one,                          &
                     g_ab_dc,                      &
                     length_1)
!
         call wf%mem%dealloc(x_ib_dc, wf%n_o*length_2, length_4*length_3)
!
!        g_vv_vv = g_ab_cd(ab, cd) += g_ab_dc(ab, dc)
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
         call wf%mem%dealloc(g_ab_dc, length_1*length_2, length_4*length_3)
!
      endif
!
   end subroutine t1_transform_vv_vv_ccs
!
!
   module subroutine store_vv_vv_electronic_repulsion_ccs(wf)
!!
!!    Store vvvv Electronic Repulsion  
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
      integer(i15) :: a = 0, b = 0, c = 0, d = 0, bcd = 0, cd_packed = 0, I = 0, bcd_nonpacked = 0
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
!        Prints
!
         if (wf%settings%print_level == 'developer') write(unit_output,'(/t3,a36)')'Integrals g_abcd are stored on disk.'
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
            call wf%mem%alloc(L_ab_J, (wf%n_v)*b_length, wf%n_J)
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
               call wf%mem%alloc(L_cd_J, (wf%n_v)*d_length, wf%n_J)
!
               call wf%read_cholesky_ab(L_cd_J, 1, wf%n_v, d_first, d_last)
!
               call wf%mem%alloc(g_a_bcd, (wf%n_v), b_length*(wf%n_v)*d_length)
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
               call wf%mem%dealloc(L_cd_J, (wf%n_v)*d_length, wf%n_J)
!
!              Save the integrals to disk 
!
               do d = 1, d_length
                  do c = 1, wf%n_v
                     do b = 1, b_length
!
!                       Calculate record number 
!
                        cd_packed = index_packed(c, d + d_first - 1)
                        bcd = index_two(b + b_first - 1, cd_packed, wf%n_v) ! Packed index!
                        bcd_nonpacked = index_three(b, c, d, b_length, wf%n_v) ! Nonpacked index! 
!
                        if (c .ge. (d + d_first - 1)) then 
!
!                          Write integrals to that record 
!
                           write(unit_g_abcd, rec=bcd) (g_a_bcd(I, bcd_nonpacked), I = 1, wf%n_v)
!
                        endif
!
                     enddo
                  enddo
               enddo
!
               call wf%mem%dealloc(g_a_bcd, (wf%n_v), b_length*(wf%n_v)*d_length)
!
            enddo ! End of batches over d
!
            call wf%mem%dealloc(L_ab_J, (wf%n_v)*b_length, wf%n_J)
!
         enddo ! End of batches over b 
!
!        Test for file handling error 
!
         if (ioerror .ne. 0) then 
!
            write(unit_output,'(t3,a)') 'Error: write error in store_electronic_repulsion_integrals_ccs'
            stop
!
         endif
!
!        Close file 
!
         close(unit_g_abcd)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then
! 
            write(unit_output,'(t3,a36,f14.8)') 'Time to store g_abcd (seconds):     ', end_timer - begin_timer
            flush(unit_output)
!
         endif
!
      endif 
!
   end subroutine store_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine store_t1_vv_vv_electronic_repulsion_ccs(wf)
!!
!!    Store t1 vvvv Electronic Repulsion Integrals 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Tests whether it is possible to store t1-transformed vir-vir-vir-vir integrals and,
!!    if possible, writes the integrals to disk 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      integer  :: required_space 
      real(dp) :: required_space_gb
!
      integer(i15) :: unit_g_abcd = -1    ! g_abcd, non-transformed 
      integer(i15) :: unit_g_t1_abcd = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
      integer(i15) :: rec_number = -1     ! The record where g_abcd is positioned 
!
!     Batching variables
!
      integer(i15) :: required_mem = -1
      integer(i15) :: available_mem = -1 
!
      integer(i15) :: b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0, b_n_batch = 0
      integer(i15) :: d_first = 0, d_last = 0, d_length = 0, d_max_length = 0, d_batch = 0
!
      integer(i15) :: a = 0, b = 0, c = 0, d = 0, bcd = 0, cd_packed = 0, I = 0, bcd_nonpacked = 0
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
      required_space = (wf%n_v)**4
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
!        Open and delete file containing non-transformed integrals - if it exists, of course
!
         call generate_unit_identifier(unit_g_abcd)
         open(unit=unit_g_abcd, file='g_abcd', action='read', status='unknown', &
               access='direct', form='unformatted', recl=dp*(wf%n_v), iostat=ioerror) 
         close(unit=unit_g_abcd, status='delete')
!
         call cpu_time(begin_timer)
!
!        Open file for writing integrals - one record: (a, bcd) = (1:n_v, bcd)
!
         call generate_unit_identifier(unit_g_t1_abcd)
         open(unit=unit_g_t1_abcd, file='g_t1_abcd', action='write', status='unknown', &
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
            call wf%mem%alloc(L_ab_J, (wf%n_v)*b_length, wf%n_J)
            call wf%get_cholesky_ab(L_ab_J, 1, wf%n_v, b_first, b_last)
!
            do d_batch = 1, b_batch ! Batch over d index; restricted by b index 
!
               call batch_limits(d_first, d_last, d_batch, d_max_length, wf%n_v)
               d_length = d_last - d_first + 1 
!
!              Calculate the integrals g_ab_cd, where b and d are restricted by
!              the current batching limits 
!
               call wf%mem%alloc(L_cd_J, (wf%n_v)*d_length, wf%n_J)
!
               call wf%get_cholesky_ab(L_cd_J, 1, wf%n_v, d_first, d_last)
!
               call wf%mem%alloc(g_a_bcd, (wf%n_v), b_length*(wf%n_v)*d_length)
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
               call wf%mem%dealloc(L_cd_J, (wf%n_v)*d_length, wf%n_J)
!
!              Save the integrals to disk 
!
               do d = 1, d_length
                  do c = 1, wf%n_v
                     do b = 1, b_length
!
!                       Calculate record number 
!
                        bcd_nonpacked = index_three(b, c, d, b_length, wf%n_v) ! Nonpacked index 
!
!                       Write integrals to that record 
!
                        write(unit_g_t1_abcd, rec=bcd_nonpacked) (g_a_bcd(I, bcd_nonpacked), I = 1, wf%n_v)
!
                     enddo
                  enddo
               enddo
!
               call wf%mem%dealloc(g_a_bcd, (wf%n_v), b_length*(wf%n_v)*d_length)
!
            enddo ! End of batches over d
!
            call wf%mem%dealloc(L_ab_J, (wf%n_v)*b_length, wf%n_J)
!
         enddo ! End of batches over b 
!
!        Test for file handling error 
!
         if (ioerror .ne. 0) write(unit_output,'(t3,a)') 'Error: write error in store_t1_vv_vv_electronic_repulsion_integrals_ccs'
!
!        Close file 
!
         close(unit_g_t1_abcd)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then
! 
            write(unit_output,'(t6,a39,f14.8)') 'Time used to store t1 g_abcd (seconds):', end_timer - begin_timer
            flush(unit_output)
!
         endif
!
      endif 
!
   end subroutine store_t1_vv_vv_electronic_repulsion_ccs
!
!
   module subroutine store_t1_vo_ov_electronic_repulsion_ccs(wf)
!!
!!    Store t1 voov Electronic Repulsion (CCS) 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Tests whether it is possible to store t1-transformed vir-occ-occ-vir integrals and,
!!    if possible, writes the integrals to disk 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      integer  :: required_space 
      real(dp) :: required_space_gb
!
      integer(i15) :: unit_g_t1_aijb = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
!
      real(dp) :: begin_timer, end_timer
!
!     Integral
!
      real(dp), dimension(:,:), allocatable :: g_ai_jb
      real(dp), dimension(:,:), allocatable :: L_jb_J
      real(dp), dimension(:,:), allocatable :: L_ai_J
!
      character(len=40) :: integral_type 
!
      integer(i15) :: a = 0, i = 0, ai = 0
      integer(i15) :: jb = 0
!
!     Disk space required to store g_vo_ov 
!
      required_space = ((wf%n_v)**2)*((wf%n_o)**2)
!
!     This is the required space in number of double precision numbers (8 bytes per such number).
!     We convert this number to gigabytes. 
!
      required_space = 8*required_space ! in bytes
!
!     Required in giga bytes
!
      required_space_gb = real(required_space)*(1.0D-9)
!
!     Test whether there is room for the integrals & save if this is the case 
!
      if (required_space_gb .lt. wf%settings%disk_space) then 
!
         call cpu_time(begin_timer)
!
         call wf%mem%alloc(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!        Alllocate Cholesky vectors
!
         call wf%mem%alloc(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
         call wf%mem%alloc(L_jb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ai(L_ai_J)
         call wf%get_cholesky_ia(L_jb_J)
!
!        Construct integral
!
         call dgemm('N', 'T',             &
                     (wf%n_o)*(wf%n_v),   &
                     (wf%n_o)*(wf%n_v),   &
                     wf%n_J,              &
                     one,                 &
                     L_ai_J,              &
                     (wf%n_o)*(wf%n_v),   &
                     L_jb_J,              &
                     (wf%n_o)*(wf%n_v),   &
                     zero,                &
                     g_ai_jb,             &
                     (wf%n_o)*(wf%n_v))
!
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
         call wf%mem%dealloc(L_jb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!  
!        Open file for writing integrals - record (ai), record length (n_o)*(n_v) (bj)
!
         call generate_unit_identifier(unit_g_t1_aijb)
         open(unit=unit_g_t1_aijb, file='g_t1_aijb', action='write', status='unknown', &
               access='direct', form='unformatted', recl=dp*((wf%n_v)*(wf%n_o)), iostat=ioerror)
!
         if (ioerror .ne. 0) write(unit_output,*) &
         'Error: error while opening file in store_t1_vo_ov_electronic_repulsion_ccs', ioerror
!
         do a = 1, wf%n_v
            do i = 1, wf%n_o
!
!                 Calculate record number 
!
                  ai = index_two(a, i, wf%n_v)
!                  
!                 Write integrals to that record 
!
                  write(unit_g_t1_aijb, rec=ai, iostat=ioerror) (g_ai_jb(ai, jb), jb = 1, wf%n_v*wf%n_o)
!
            enddo
         enddo
!
         if (ioerror .ne. 0) write(unit_output,'(t3,a)') &
            'Error: write error in store_t1_vo_ov_electronic_repulsion_integrals_ccs'
!
         call wf%mem%dealloc(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then
! 
            write(unit_output,'(t6,a39,f14.8)') 'Time used to store t1 g_aijb (seconds):', end_timer - begin_timer
            flush(unit_output)
!
         endif
         close(unit_g_t1_aijb)
!
      endif
!
   end subroutine store_t1_vo_ov_electronic_repulsion_ccs
!
!
   module subroutine store_t1_vv_vo_electronic_repulsion_ccs(wf)
!!
!!    Store t1 vvvo Electronic Repulsion (CCS) 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Tests whether it is possible to store t1-transformed vir-vir-vir-occ integrals and,
!!    if possible, writes the integrals to disk 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      integer  :: required_space 
      real(dp) :: required_space_gb
!
      integer(i15) :: unit_g_t1_abci = -1 ! g_abic, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
!
      integer(i15) :: required_mem = 0, available_mem = 0, a_batch = 0, bci = 0, a_full = 0, a = 0, b = 0, ci = 0
      integer(i15) :: a_n_batch = 0, a_max_length = 0, a_length = 0, a_first = 0, a_last = 0, ab_rec = 0, ab = 0
!
      real(dp) :: begin_timer, end_timer
!
!     Integral
!
      real(dp), dimension(:,:), allocatable :: g_ab_ci
!
      character(len=40) :: integral_type 
!
      integer(i15) :: a = 0, i = 0, ai = 0
      integer(i15) :: jb = 0
!
!     Disk space required to store g_vo_vv 
!
      required_space = ((wf%n_v)**3)*((wf%n_o)**2)
!
!     This is the required space in number of double precision numbers (8 bytes per such number).
!     We convert this number to gigabytes. 
!
      required_space = 8*required_space ! in bytes
!
!     Required in giga bytes
!
      required_space_gb = real(required_space)*(1.0D-9)
!
!     Test whether there is room for the integrals & save if this is the case 
!
      if (required_space_gb .lt. wf%settings%disk_space) then 
!
         call cpu_time(begin_timer)
!  
!        Open file for writing integrals
!
         call generate_unit_identifier(unit_g_t1_abci)
         open(unit=unit_g_t1_abci, file='g_t1_abci', action='write', status='unknown', &
               access='direct', form='unformatted', recl=dp*((wf%n_o)*(wf%n_v)), iostat=ioerror)
!
         if (ioerror .ne. 0) write(unit_output,*) &
         'Error: error while opening file in store_t1_vv_vo_electronic_repulsion_ccs', ioerror
!
!        In calculating g_ab_ic, we will batch over the a index 
!
         required_mem = max((wf%n_v)**2*(wf%n_J) + 2*(wf%n_v)*(wf%n_o)*(wf%n_J), & ! Needed to get L_ab_J
                           ((wf%n_v)**3)*(wf%n_o) + 2*(wf%n_v)*(wf%n_o)*(wf%n_J)) ! Needed to get g_ab_ci
!         
         required_mem  = 4*required_mem
         available_mem = get_available()
!
         a_max_length = 0
         call num_two_batch(required_mem, available_mem, a_max_length, a_n_batch, wf%n_v)
!
         do a_batch = 1, a_n_batch ! Batch over a index 
!
            call batch_limits(a_first, a_last, a_batch, a_max_length, wf%n_v)
            a_length = a_last - a_first + 1  

            call wf%mem%alloc(g_ab_ci, a_length*(wf%n_v), (wf%n_o)*(wf%n_v))
!
            call wf%get_vv_vo_electronic_repulsion(g_ab_ci, a_first, a_last, 1, wf%n_v, 1, wf%n_v, 1, wf%n_o)
!
            do a = 1, a_length
               do b = 1, wf%n_v
!
!                 Calculate record number 
!
                  a_full = a + a_first - 1
                  ab_rec = index_two(a_full, b, wf%n_v)
!
!                 Calculate ab index for (possibly restricted) g_ab_ci integral
!
                  ab = index_two(a, b, a_length)                  
!
!                 Write integrals to that record 
!
                  write(unit_g_t1_abci, rec=ab_rec, iostat=ioerror) (g_ab_ci(ab, ci), ci = 1, (wf%n_o)*(wf%n_v))
!
               enddo
            enddo
!
            if (ioerror .ne. 0) write(unit_output,'(t3,a)') &
               'Error: write error in store_t1_vv_vo_electronic_repulsion_integrals_ccs'
!
            call wf%mem%dealloc(g_ab_ci, a_length*(wf%n_v), (wf%n_o)*(wf%n_v))
!
         enddo ! End of batches over a 
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then
! 
            write(unit_output,'(t6,a39,f14.8)') 'Time used to store t1 g_abci (seconds):', end_timer - begin_timer
            flush(unit_output)
!
         endif

         close(unit_g_t1_abci)
!
      endif
!
   end subroutine store_t1_vv_vo_electronic_repulsion_ccs
!
!
   module subroutine read_t1_vo_ov_electronic_repulsion_ccs(wf, x_vo_ov, & 
                                             index1_first, index1_last,  &
                                             index2_first, index2_last,  &
                                             index3_first, index3_last,  &
                                             index4_first, index4_last)
!!
!!    Read t1 voov Electronic Repulsion Integrals 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Reads the T1-transformed vir-occ-occ-vir integrals from file.
!!
!!    g_t1_ai_jb is written such that we have full flexibility with respect to ai indices. However, get_vo_ov_ccs has
!!    full flexibility wrt. all indices and we must check wether b and j are full space or not
!!    (they will presumably always be full space indices).
!!
      implicit none 
!
      class(ccs) :: wf
!
!     Integral
!
      real(dp), dimension(:,:) :: x_vo_ov 
!
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last

!
      integer(i15) :: unit_g_t1_aijb = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
!
      real(dp) :: begin_timer, end_timer
!
      real(dp), dimension(:,:), allocatable :: g_ai_jb_full 
!
      character(len=40) :: integral_type 
!
      integer(i15) :: a = 0, i = 0, ai = 0
      integer(i15) :: b = 0, j = 0, jb = 0, jb_full = 0
!
      integer(i15) :: length_a = 0, length_i = 0, length_j = 0, length_b = 0
!
!     Lengths
!
      length_a = index1_last - index1_first + 1
      length_i = index2_last - index2_first + 1
      length_j = index3_last - index3_first + 1
      length_b = index4_last - index4_first + 1
!

      call cpu_time(begin_timer)
!  
!     Open file for reading integrals - record (ai), record length (n_o)*(n_v) (bj)
!
      call generate_unit_identifier(unit_g_t1_aijb)
      open(unit=unit_g_t1_aijb, file='g_t1_aijb', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*((wf%n_v)*(wf%n_o)), iostat=ioerror)
!
      if (ioerror .ne. 0) write(unit_output,'(t3,a)') &
      'Error: error while opening file in read_t1_vo_ov_electronic_repulsion_ccs'
!  
!     Check if length_b = wf%n_v and length_j = wf%n_o -> this will usually be the case
!
      if (length_b .eq. wf%n_v .and. length_j .eq. wf%n_o) then
!
         do a = 1, length_a
            do i = 1, length_i
!
!              Calculate record number 
!
               ai = index_two(a + index1_first - 1, i + index2_first - 1, wf%n_v)
!               
!              Write integrals to that record 
!
               read(unit_g_t1_aijb, rec=ai, iostat=ioerror) (x_vo_ov(ai, jb), jb = 1, (wf%n_v)*(wf%n_o))
!
            enddo
         enddo
!
      else ! -> Unlikely to happen
!
!        Must limit jb index, will first read for all jb's
!
         call wf%mem%alloc(g_ai_jb_full, length_a*length_i, (wf%n_v)*(wf%n_o))
!
         do a = 1, length_a
            do i = 1, length_i
!
!              Calculate record number 
!
               ai = index_two(a + index1_first - 1, i + index2_first - 1, wf%n_v)
!               
!              Write integrals to that record 
!
               read(unit_g_t1_aijb, rec=ai, iostat=ioerror) (g_ai_jb_full(ai, jb), jb = 1, (wf%n_v)*(wf%n_o))
!
            enddo
         enddo
!
!        Then sort away jb elements not required - pretty sure we will never be used
!
         do j = 1, length_j
            do b = 1, length_b
!
               jb_full = index_two(j + index3_first - 1, b + index4_first - 1, wf%n_o)
               jb = index_two(j, b, length_j)
!
               x_vo_ov(:, jb) = g_ai_jb_full(:, jb_full)
!
            enddo
         enddo
!
         call wf%mem%dealloc(g_ai_jb_full, length_a*length_i, (wf%n_v)*(wf%n_o))
!
      endif
!
      if (ioerror .ne. 0) write(unit_output,'(t3,a)') &
         'Error: read error in read_t1_vo_ov_electronic_repulsion_integrals_ccs'
      close(unit_g_t1_aijb)
!
      call cpu_time(end_timer)
!
      if (wf%settings%print_level == 'developer') then
! 
         ! write(unit_output,'(t3,a39,f14.8)') 'Time used to read t1 g_aijb (seconds):', end_timer - begin_timer
         ! flush(unit_output)
!
      endif
!
   end subroutine read_t1_vo_ov_electronic_repulsion_ccs
!
!
   module subroutine store_t1_vv_ov_electronic_repulsion_ccs(wf)
!!
!!    Store t1 ovvv Electronic Repulsion Integrals 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Tests whether it is possible to store t1-transformed occ-vir-vir-vir integrals and,
!!    if possible, writes the integrals to disk 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      integer  :: required_space 
      real(dp) :: required_space_gb
!
      integer(i15) :: unit_g_t1_bcia = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
      integer(i15) :: rec_number = -1     ! The record where g_abcd is positioned 
!
      real(dp) :: begin_timer, end_timer
!
!     Batching variables
!
      integer(i15) :: required_mem = -1
      integer(i15) :: available_mem = -1 
!
      integer(i15) :: b_first = 0, b_last = 0, b_length = 0, b_max_length = 0, b_batch = 0, b_n_batch = 0
!
      integer(i15) :: a = 0, b = 0, c = 0, bc = 0, ia = 0, record_number = 0
!
!     The electronic repulsion integral 
!
      real(dp), dimension(:,:), allocatable :: g_bc_ia ! g_iabc 
      real(dp), dimension(:,:), allocatable :: L_bc_J
      real(dp), dimension(:,:), allocatable :: L_ia_J
!
!     Calculate the disk space (in GB) required to store the occ-vir-vir-vir integrals
!
      required_space = ((wf%n_v)**3)*(wf%n_o)
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
!        Begin timings
!
         call cpu_time(begin_timer)
!
!        Open file for writing integrals - one record
!
         call generate_unit_identifier(unit_g_t1_bcia)
         open(unit=unit_g_t1_bcia, file='g_t1_bcia', action='write', status='unknown', &
               access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
!        In calculating g_ab_cd, we will batch over the b and d indices 
!
         required_mem = max(2*(wf%n_v)**2*(wf%n_J) + 4*(wf%n_v)*(wf%n_o)*(wf%n_J), & ! Needed to get L_ia_J & L_bc_J
                           (wf%n_v)**3*(wf%n_o) + (wf%n_v)**2*(wf%n_J) + (wf%n_v)*(wf%n_o)*(wf%n_J)) ! Needed to get g_ia_bc
!         
         required_mem  = 4*required_mem
         available_mem = get_available()
!
         b_max_length = 0
         call num_batch(required_mem, available_mem, b_max_length, b_n_batch, wf%n_v)
!
         do b_batch = 1, b_n_batch ! Batch over b index 
!
            call batch_limits(b_first, b_last, b_batch, b_max_length, wf%n_v)
            b_length = b_last - b_first + 1  
!
            call wf%mem%alloc(g_bc_ia, b_length*(wf%n_v), (wf%n_v)*(wf%n_o))
!
            call wf%mem%alloc(L_bc_J, b_length*(wf%n_v), wf%n_J)
            call wf%mem%alloc(L_ia_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!        Get T1-transformed Cholesky vectors
!
         call wf%get_cholesky_ab(L_bc_J,b_first, b_last, 1, wf%n_v)
         call wf%get_cholesky_ia(L_ia_J)
!
!        Construct integral
!
         call dgemm('N', 'T',             &
                     b_length*(wf%n_v),   &
                     (wf%n_v)*(wf%n_o),   &
                     wf%n_J,              &
                     one,                 &
                     L_bc_J,              &
                     b_length*(wf%n_v),   &
                     L_ia_J,              &
                     (wf%n_v)*(wf%n_o),   &
                     zero,                &
                     g_bc_ia,             &
                     b_length*(wf%n_v))
!
!        Deallocate Cholesky vectors
!
         call wf%mem%dealloc(L_bc_J, b_length*(wf%n_v), wf%n_J)
         call wf%mem%dealloc(L_ia_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
!           Save the integrals to disk 
!
            do b = b_first, b_last
               do c = 1, wf%n_v
!
!                 Write integrals to that record 
!
                  bc = index_two(b - b_first + 1, c, b_length)
                  record_number = index_two(b, c, wf%n_v) 
                  write(unit_g_t1_bcia, rec=record_number) (g_bc_ia(bc, ia), ia = 1, (wf%n_v)*(wf%n_o))
!
               enddo
            enddo
!
!
         enddo ! End of batches over b 
!
!        Test for file handling error 
!
         if (ioerror .ne. 0) write(unit_output,'(t3,a)') 'Error: write error in store_t1_vv_ov_electronic_repulsion_integrals_ccs'
!
!        Close file 
!
         close(unit_g_t1_bcia)
!
         call cpu_time(end_timer)
!
         if (wf%settings%print_level == 'developer') then
! 
            write(unit_output,'(t6,a39,f14.8)') 'Time used to store t1 g_bcia (seconds):', end_timer - begin_timer
            flush(unit_output)
!
         endif
!
      endif 
!
   end subroutine store_t1_vv_ov_electronic_repulsion_ccs
!
!
   module subroutine read_t1_vv_ov_electronic_repulsion_ccs(wf, x_vv_ov,& 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Read t1 vvov Electronic Repulsion Integrals 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Tests whether it is possible to store t1-transformed vir-vir-occ-vir integrals and,
!!    if possible, writes the integrals to disk 
!!
!!    Assumes batching over either b or c, so please don't batch over a because then routine will not work.
!!
!!
      implicit none 
!
      class(ccs) :: wf
!
!     Integral
!
      real(dp), dimension(:,:) :: x_vv_ov
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
      real(dp) :: begin_timer, end_timer
!
      integer(i15) :: unit_g_t1_bcia = -1 ! g_abcd, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
      integer(i15) :: rec_number = -1     ! The record where g_abcd is positioned 
!
      integer(i15) :: length_a = 0, length_i = 0, length_c = 0, length_b = 0
!
      integer(i15) :: a = 0, b = 0, c = 0, bc = 0, record_number = 0, i = 0, ia = 0, ia_full = 0
!
      real(dp), dimension(:,:), allocatable :: g_bc_ia 
!
!     Lengths
!
      length_b = index1_last - index1_first + 1
      length_c = index2_last - index2_first + 1
      length_i = index3_last - index3_first + 1
      length_a = index4_last - index4_first + 1
!
!     Begin timings
!
      call cpu_time(begin_timer)
!
!     Open file for writing integrals - one record
!
      call generate_unit_identifier(unit_g_t1_bcia)
      open(unit=unit_g_t1_bcia, file='g_t1_bcia', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
      do b = index1_first, index1_last
         do c = index2_first, index2_last
!
!           Write integrals to that record 
!
            bc = index_two(b - index1_first + 1, c - index2_first + 1, length_b)
            record_number = index_two(b, c, wf%n_v) 
            read(unit_g_t1_bcia, rec=record_number) (x_vv_ov(bc, ia), ia = 1, (wf%n_v)*(wf%n_o))

!
         enddo
      enddo
!
!     Test for file handling error 
!
      if (ioerror .ne. 0) write(unit_output,'(t3,a)') 'Error: write error in adre_t1_vv_ov_electronic_repulsion_integrals_ccs'
!
!     Close file 
!
      close(unit_g_t1_bcia)
!
   end subroutine read_t1_vv_ov_electronic_repulsion_ccs
!
!
   module subroutine read_t1_vv_vo_electronic_repulsion_ccs(wf, x_vv_vo,& 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!    Read t1 vvvo Electronic Repulsion Integrals 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Reads the T1-transformed vir-vir-vir-occ integrals from file.
!!
!!    The integrals are stored on file as (ab,ci) = (ab, :), where a 
!!    is the record number and : denotes all the bci elements.
!!
!!    The recommended use is therefore to batch over the a or b index,
!!    as this will involve the no wasteful read statements
!!
      implicit none 
!
      class(ccs) :: wf 
!
!     Integral
!
      real(dp), dimension(:,:) :: x_vv_vo
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
      integer(i15) :: unit_g_t1_abci = -1 ! g_abci, electronic repulsion integrals  
      integer(i15) :: ioerror = -1        ! Error integer for file handling
!
      real(dp) :: begin_timer, end_timer
!
      real(dp), dimension(:,:), allocatable :: g_AB_ci ! Holds the elements AB, ci for a given AB index 
!
      character(len=40) :: integral_type 
!
      integer(i15) :: a = 0, a_full = 0, b_full = 0, ab_rec = 0, c = 0, ci = 0, i = 0, ci_full = 0, b = 0, ab = 0
!
      integer(i15) :: length_a = 0, length_b = 0, length_c = 0, length_i = 0
!
!     Lengths
!
      length_a = index1_last - index1_first + 1
      length_b = index2_last - index2_first + 1
      length_c = index3_last - index3_first + 1
      length_i = index4_last - index4_first + 1
!
      call cpu_time(begin_timer)
!  
!     Open file for reading integrals - record (ai), record length (n_o)*(n_v) (bj)
!
      call generate_unit_identifier(unit_g_t1_abci)
      open(unit=unit_g_t1_abci, file='g_t1_abci', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*((wf%n_o)*(wf%n_v)), iostat=ioerror)
!
      if (ioerror .ne. 0) write(unit_output,'(t3,a)') &
      'Error: error while opening file in read_t1_vv_vo_electronic_repulsion_ccs'
!
      call wf%mem%alloc(g_AB_ci, 1, (wf%n_o)*(wf%n_v))
      g_AB_ci = zero
!
      do a = 1, length_a
         do b = 1, length_b
!
!           For a given A and B, read into g_AB_ci 
!
            a_full = a + index1_first - 1
            b_full = b + index2_first - 1
!
            ab_rec = index_two(a_full, b_full, wf%n_v)
!
            read(unit_g_t1_abci, rec=ab_rec, iostat=ioerror) (g_AB_ci(1,ci_full), ci_full = 1, (wf%n_o)*(wf%n_v))
!
!           Place the result into the incoming integral array (ab_ci)
!
            do i = 1, length_i
               do c = 1, length_c
!
                  ci = index_two(c, i, length_c)
!
                  ab = index_two(a, b, length_a)
!
                  ci_full = index_two(c + index3_first - 1,  &
                                       i + index4_first - 1, &
                                       wf%n_v)
!
                  x_vv_vo(ab, ci) = g_AB_ci(1, ci_full)
!
               enddo
            enddo
         enddo
!
      enddo ! End of read loop over a 
!
      call wf%mem%dealloc(g_AB_ci, 1, (wf%n_o)*(wf%n_v))
!
      if (ioerror .ne. 0) write(unit_output,'(t3,a)') &
         'Error: read error in read_t1_vv_vo_electronic_repulsion_integrals_ccs'
!
      close(unit_g_t1_abci)
!
      call cpu_time(end_timer)
!
      if (wf%settings%print_level == 'developer') then
! 
         ! write(unit_output,'(t3,a39,f14.8)') 'Time used to read t1 g_abci (seconds):', end_timer - begin_timer
         ! flush(unit_output)
!
      endif
!
   end subroutine read_t1_vv_vo_electronic_repulsion_ccs
!
!
   module function get_vvvv_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!!
!!    Get vvvv required memory (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
!!    Calculates and returns required memory to make vvvv electronic repulsion integral.
!!
!!    dim_1, dim_2, dim_3, and dim_4 are the full dimension of index 1-4.
!!    They will typically be wf%n_v and are therefore optionals, however will not be wf%n_v for ML 
!!
      implicit none
!
      class(ccs), intent(in)              :: wf 
!  
      integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
      integer(i15) :: get_vvvv_required_mem_ccs
!
      logical :: vvvv_on_file    = .false.
      logical :: t1_vvvv_on_file = .false.
!
      if (present(dim_1) .and. present(dim_2) .and. present(dim_3) .and. present(dim_4)) then
!
         get_vvvv_required_mem_ccs = (dim_1*dim_2*dim_3*dim_4)
!
!        Check if vvvv integrals are on file
!
         inquire(file='g_abcd',exist=vvvv_on_file)
         inquire(file='g_t1_abcd',exist=t1_vvvv_on_file)
!
         if ( t1_vvvv_on_file) then
!
!           We are reading T1 transformed vvvv electronic repulsion integral, 
!           this only requires the size of the array itself
!
            get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs*dp
            return
!
         elseif (vvvv_on_file) then
!
!           We are reading vvvv electronic repulsion integral, 
!           this only requires the size of the array itself plus T1 transformation
!
            get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs + &
                     max((dim_2)*(wf%n_o)*(wf%n_J) + (dim_3)*(dim_4)*(wf%n_J) + (dim_2)*(dim_3)*(dim_4)*(wf%n_o)       , &
                         (dim_4)*(wf%n_o)*(wf%n_J) + (dim_1)*(dim_2)*(dim_4)*(wf%n_o) + (dim_1)*(dim_2)*(wf%n_J)       , &
                         (dim_4)*(wf%n_o)*(wf%n_J) + (dim_1)*(dim_2)*(dim_4)*(wf%n_o) + (dim_1)*(dim_2)*(dim_3)*(dim_4), &
                         (dim_4)*(wf%n_o)*(wf%n_J) + (dim_2)*(dim_4)*(wf%n_o**2) &
                              + (dim_2)*(wf%n_o)*(wf%n_J) + (dim_2)*(dim_4)*(dim_3)*(dim_4), &
                         (dim_2)*(dim_4)*(wf%n_o**2) + (dim_2)*(dim_3)*(dim_4)*(wf%n_o) + (dim_2)*(dim_4)*(dim_3)*(dim_4))
!
            get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs*dp
            return
!
         endif
!
      elseif (.not. (present(dim_1) .and. present(dim_2) .and. present(dim_3) .and. present(dim_4))) then
!
         get_vvvv_required_mem_ccs = (wf%n_v**4)*dp
!
!        Check if vvvv integrals are on file
!
         inquire(file='g_abcd',exist=vvvv_on_file)
         inquire(file='g_t1_abcd',exist=t1_vvvv_on_file)
!
         if ( t1_vvvv_on_file) then
!
!           We are reading T1 transformed vvvv electronic repulsion integral, 
!           this only requires the size of the array itself
!
            get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs*dp
            return
!
         elseif (vvvv_on_file) then
!
!           We are reading vvvv electronic repulsion integral, 
!           this only requires the size of the array itself plus T1 transformation
!
            get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs + &
                      max((wf%n_v)*(wf%n_o)*(wf%n_J) + (wf%n_v**2)*(wf%n_J) + (wf%n_v**3)*(wf%n_o), &
                          (wf%n_v)*(wf%n_o)*(wf%n_J) + (wf%n_v**4) + (wf%n_v**3)*(wf%n_o)         , &
                          2*(wf%n_v)*(wf%n_o)*(wf%n_J) + (wf%n_v**2)*(wf%n_o**2)                  , &
                          (wf%n_v**2)*(wf%n_o**2) + (wf%n_v**3)*(wf%n_o)                          , &
                          (wf%n_v**4) + (wf%n_v**3)*(wf%n_o)) 
!
            get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs*dp
            return
!
         endif
!
      else
!
         write(unit_output,*) 'Error: call to get_vvvv_required_mem is missing some arguments'
         stop
!
      endif
!
!     We are constructing the integral from T1-transformed Cholesky vectors
!
      get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs + (dim_1)*(dim_2)*(wf%n_J)
!
      get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs + max(2*(wf%n_o)*(dim_2)*(wf%n_J), &
                                       (wf%n_o)*(dim_2)*(wf%n_J) + (dim_1)*(dim_2)*(wf%n_J), &
                                       (dim_3)*(dim_4)*(wf%n_J) + 2*(wf%n_o)*(dim_4)*(wf%n_J), &
                                       2*(dim_3)*(dim_4)*(wf%n_J) + (wf%n_o)*(dim_4)*(wf%n_J))
!
      get_vvvv_required_mem_ccs = get_vvvv_required_mem_ccs*dp
!
   end function get_vvvv_required_mem_ccs
!
!
   integer function get_vvvo_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!!
!!    Get vvvo required memory (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
!!    Calculates and returns required memory to make vvvv electronic repulsion integral.
!!
!!    dim_1, dim_2, dim_3, and dim_4 are the full dimension of index 1-4.
!!    They will typically be wf%n_v for 1-3 and wf%n_o for 4 and are therefore optionals,  
!!    however will not necessarily have these values for ML 
!!
      implicit none
!
      class(ccs), intent(in)              :: wf 
!  
      integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
      integer(i15) :: get_vvvo_required_mem_ccs
!
   end function get_vvvo_required_mem_ccs
!
!
   integer function get_vvov_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!!
!!    Get vvvo required memory (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
!!    Calculates and returns required memory to make vvvv electronic repulsion integral.
!!
!!    dim_1, dim_2, dim_3, and dim_4 are the full dimension of index 1-4.
!!    They will typically be wf%n_v for 1, 2, and 4 and wf%n_o for 3 and are therefore optionals,  
!!    however will not necessarily have these values for ML 
!!
      implicit none
!
      class(ccs), intent(in)              :: wf 
!  
      integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
      integer(i15) :: get_vvvo_required_mem_ccs
!
   end function get_vvov_required_mem_ccs
!
!
   integer function get_vvoo_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!!
!!    Get vvvo required memory (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
!!    Calculates and returns required memory to make vvvv electronic repulsion integral.
!!
!!    dim_1, dim_2, dim_3, and dim_4 are the full dimension of index 1-4.
!!    They will typically be wf%n_v for 1 and 2 and wf%n_o for 3 and 4 and are therefore optionals,  
!!    however will not necessarily have these values for ML 
!!
      implicit none
!
      class(ccs), intent(in)              :: wf 
!  
      integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
      integer(i15) :: get_vvvo_required_mem_ccs
!
   end function get_vvoo_required_mem_ccs
!
!
end submodule integrals