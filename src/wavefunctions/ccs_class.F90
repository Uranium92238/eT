module ccs_class
!
!!
!!    Coupled cluster singles (ccs) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
   use hf_class
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
   subroutine initialize_ccs(wf, ref_wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      class(hf) :: ref_wf
!
      wf%name = 'ccs'
      write(output%unit, *)'init'
      flush(output%unit)
!
      wf%system = ref_wf%system
!
      wf%n_ao   = ref_wf%n_ao
      wf%n_mo   = ref_wf%n_mo
      wf%n_o    = ref_wf%n_o
      wf%n_v    = ref_wf%n_v
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
