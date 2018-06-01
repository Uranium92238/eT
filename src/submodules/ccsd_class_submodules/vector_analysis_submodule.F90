submodule (ccsd_class) vector_analysis
!!
!!    Vector analysis submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2018
!!
!!    Contains routines that analyze vectors, which includes general analysis
!!    of, e.g., the ground and excited state vectors.
!!
contains
!
!
   module subroutine print_dominant_doubles_ccsd(wf, vec)
!!
!!    Print dominant singles (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:, :) :: vec ! Assumed to have n_t2am number of elements
!
      real(dp), dimension(:,:), allocatable :: sorted_max_vec_doubles
      integer(i15), dimension(:,:), allocatable :: index_list_doubles
!
      integer(i15) :: i = 0
!
      call wf%mem%alloc(sorted_max_vec_doubles, min(wf%n_t2am, 20), 1)
      call wf%mem%alloc_int(index_list_doubles, min(wf%n_t2am, 20), 4)
!
      call wf%analyze_double_excitation_vector(vec, min(wf%n_t2am, 20), &
                                               sorted_max_vec_doubles, index_list_doubles)
!
      write(unit_output,'(/t6,a)') 'Largest contributions to doubles vector:'
!
      write(unit_output,'(t6,a54)')'------------------------------------------------------'
      write(unit_output,'(t6,a3, 8x, a3, 8x, a3, 8x, a3, 8x, a10)')'a','i','b','j', 'amplitude'
      write(unit_output,'(t6,a54)')'------------------------------------------------------'
!
      do i = 1, min(wf%n_t2am, 20)
!
         if (abs(sorted_max_vec_doubles(i, 1)) .lt. 1.0D-03) then
!
            exit ! Never print very, very small amplitudes
!
         else
!
            write(unit_output,'(t6,i3, 8x,i3, 8x,i3, 8x, i3, 10x, f8.5)')         &
                                                      index_list_doubles(i, 1),   &
                                                      index_list_doubles(i, 2),   &
                                                      index_list_doubles(i, 3),   &
                                                      index_list_doubles(i, 4),   &
                                                      sorted_max_vec_doubles(i, 1)
         endif
!
      enddo
!
      write(unit_output,'(t6,a54)')'------------------------------------------------------'
!
   end subroutine print_dominant_doubles_ccsd
!
!
end submodule vector_analysis
