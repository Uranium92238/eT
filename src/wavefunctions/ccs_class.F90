module ccs_class
!
!!
!!    Coupled cluster singles (ccs) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
!
   use reordering
   use array_utilities
   use array_analysis
   use interval_class
   use index
!
   implicit none
!
   type, extends(wavefunction):: ccs
!
   contains
!
!     Initialize and finalize wavefunction
!
      procedure :: initialize => initialize_ccs
      procedure :: finalize   => finalize_ccs
!
!
   end type ccs
!
!
contains
!
subroutine initialize_ccs(wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      wf%name = 'ccs'
!
   end subroutine initialize_ccs
!
!
   subroutine finalize_ccs(wf)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
!     Nothing here yet
!
   end subroutine finalize_ccs
!
end module ccs_class
