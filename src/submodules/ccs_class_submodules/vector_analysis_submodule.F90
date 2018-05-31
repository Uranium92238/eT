submodule (ccs_class) vector_analysis
!!
!!    Vector analysis submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2018
!!
!!    Contains routines that analyze vectors, which includes general analysis
!!    of, e.g., the ground and excited state vectors.
!!
contains
!
!
   module subroutine print_dominant_singles_ccs(wf, vec)
!!
!!    Print dominant singles (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Prints the 20 largest excitations in the vector sent to the routine.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:, :) :: vec ! Assumed to have n_o*n_v number of elements
!
      real(dp) :: ddot, norm
!
      integer(i15) :: i = 0
!
      real(dp), dimension(:,:), allocatable :: sorted_max_vec_singles
      integer(i15), dimension(:,:), allocatable :: index_list_singles
!
!     Calculate the contribution from single excitations
!
      norm = sqrt(ddot(wf%n_t1am, vec, 1, vec,1))
      write(unit_output,'(/t6,a,f6.4)')'Single excitation contribution: ', norm
!
!     Analysis of vector
!
      write(unit_output,'(/t6,a)') 'Largest contributions to singles vector:'
!
      write(unit_output,'(t6,a32)')'------------------------------------------------------'
      write(unit_output,'(t6,a3, 8x, a3, 8x, a10)')'a', 'i', 'amplitude'
      write(unit_output,'(t6,a32)')'------------------------------------------------------'
!
!     Get the 20 highest amplitudes (or, if there are fewer, all of them)
!
      call wf%mem%alloc(sorted_max_vec_singles, min(wf%n_t1am, 20), 1)
      call wf%mem%alloc_int(index_list_singles, min(wf%n_t1am, 20), 2)
!
      call wf%analyze_single_excitation_vector(vec, min(wf%n_t1am, 20), &
                                               sorted_max_vec_singles, index_list_singles)
!
!     ... and then print them
!
      do i = 1, min(wf%n_t1am, 20)
!
         if (abs(sorted_max_vec_singles(i, 1)) .lt. 1.0D-03) then
!
            exit ! Never print very, very small amplitudes
!
         else
!
            write(unit_output,'(t6,i3, 8x,i3, 10x, f8.5)') index_list_singles(i, 1), &
                                                           index_list_singles(i, 2), &
                                                           sorted_max_vec_singles(i, 1)
!
         endif
!
      enddo
!
      write(unit_output,'(t6,a32)')'------------------------------------------------------'
!
   end subroutine print_dominant_singles_ccs
!
!
end submodule vector_analysis
