!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
submodule (abstract_doubles_class) file_handling_abstract_doubles
!
!!
!!    File handling submodule
!!
!!    Gathers routines that save wavefunction parameters to file,
!!    and reads them from file, plus other routines related to the 
!!    handling of the files that belong to the wavefunction.
!!
!
   implicit none
!
!
contains
!
   module subroutine read_doubles_vector_abstract_doubles(wf, X, file_)
!!
!!    Read doubles vector X in a file
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    NB: Currently only works for t and tbar-file
!!        Because the excited states use sequential_file_array which
!!        does not support several records yet.
!!
!!    Files are written with the singles part in the first record
!!    and the doubles part in the second. Thus, skip first record
!!    and read second.
!!
      implicit none 
!
      class(abstract_doubles), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(out) :: X 
!
      type(sequential_file) :: file_
!
      call file_%open_('read', 'rewind')
!
!     Skip first record because it contains only the singles part
      call file_%skip(1)
!
      call file_%read_(X, wf%n_t2)
!
      call file_%close_()
!
   end subroutine read_doubles_vector_abstract_doubles
!
end submodule file_handling_abstract_doubles 
!