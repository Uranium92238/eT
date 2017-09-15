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
      real(dp), dimension((wf%n_o)**2, (wf%n_o)**2) :: x_oo_oo
!
      integer(i15) :: index1_first, index1_last
      integer(i15) :: index2_first, index2_last
      integer(i15) :: index3_first, index3_last
      integer(i15) :: index4_first, index4_last
!
   end subroutine get_oo_oo_ccs
!
!
end submodule integrals