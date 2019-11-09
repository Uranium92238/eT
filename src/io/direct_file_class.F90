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
module direct_file_class
!
!!
!!    Direct access file class module
!!    Written by Rolf H. Myhre, May 2019
!!
!
   use kinds    
   use abstract_other_file_class, only : abstract_other_file
   use global_out, only : output
!
   type, extends(abstract_other_file) :: direct_file
!
      integer, private  :: record_dim     ! Number of words per record
      integer, private  :: word_size      ! Size of a word, default is double precision
      integer, private  :: record_length  ! record_dim*word_size
!
   contains
!
!     Open and close
!
      procedure, public :: open_    => open_direct_file
!
!
!     Write routines
!
      procedure, public :: write_r_direct_file
      procedure, public :: write_r_1_direct_file
      procedure, public :: write_r_2_direct_file
      procedure, public :: write_r_3_direct_file
      procedure, public :: write_r_4_direct_file
!
      procedure, public :: write_c_direct_file
      procedure, public :: write_c_1_direct_file
      procedure, public :: write_c_2_direct_file
      procedure, public :: write_c_3_direct_file
      procedure, public :: write_c_4_direct_file
!
      procedure, public :: write_i_direct_file
      procedure, public :: write_i_1_direct_file
      procedure, public :: write_i_2_direct_file
      procedure, public :: write_i_3_direct_file
      procedure, public :: write_i_4_direct_file
!
      generic           :: write_ => write_r_direct_file,   &
                                     write_r_1_direct_file, &
                                     write_r_2_direct_file, &
                                     write_r_3_direct_file, &
                                     write_r_4_direct_file, &
                                     write_c_direct_file,   &
                                     write_c_1_direct_file, &
                                     write_c_2_direct_file, &
                                     write_c_3_direct_file, &
                                     write_c_4_direct_file, &
                                     write_i_direct_file,   &
                                     write_i_1_direct_file, &
                                     write_i_2_direct_file, &
                                     write_i_3_direct_file, &
                                     write_i_4_direct_file
!
!     Read routines
!
      procedure, public :: read_r_direct_file
      procedure, public :: read_r_1_direct_file
      procedure, public :: read_r_2_direct_file
      procedure, public :: read_r_3_direct_file
      procedure, public :: read_r_4_direct_file
!
      procedure, public :: read_c_direct_file
      procedure, public :: read_c_1_direct_file
      procedure, public :: read_c_2_direct_file
      procedure, public :: read_c_3_direct_file
      procedure, public :: read_c_4_direct_file
!
      procedure, public :: read_i_direct_file
      procedure, public :: read_i_1_direct_file
      procedure, public :: read_i_2_direct_file
      procedure, public :: read_i_3_direct_file
      procedure, public :: read_i_4_direct_file
      generic           :: read_ => read_r_direct_file,   &
                                    read_r_1_direct_file, &
                                    read_r_2_direct_file, &
                                    read_r_3_direct_file, &
                                    read_r_4_direct_file, &
                                    read_c_direct_file,   &
                                    read_c_1_direct_file, &
                                    read_c_2_direct_file, &
                                    read_c_3_direct_file, &
                                    read_c_4_direct_file, &
                                    read_i_direct_file,   &
                                    read_i_1_direct_file, &
                                    read_i_2_direct_file, &
                                    read_i_3_direct_file, &
                                    read_i_4_direct_file
!
      final :: destructor
!
!
      procedure :: write_chimera => write_chimera_direct_file
      procedure :: read_chimera  => read_chimera_direct_file
!
   end type direct_file
!
   interface direct_file
!
      procedure new_direct_file
!
   end interface direct_file
!
contains
!
!
   function new_direct_file(name_, rec_dim, w_size) result(the_file)
!!
!!    Direct file constructer
!!    Writen by Rolf H. Myhre, May 2019
!!
!!    rec_dim is number of words in each record
!!    w_size (optional) is the size of each word, default is double precision
!!    record length is rec_dim*w_size
!!
      implicit none
!
      type(direct_file) :: the_file
!
      character(len=*), intent(in) :: name_
      integer, intent(in) :: rec_dim
      integer, intent(in), optional :: w_size
!
      if (present(w_size)) then
         if (w_size .gt. 0) then
            the_file%word_size = w_size
         else
            call output%error_msg("Word size less than zero for file "//trim(name_))
         endif
      else
         the_file%word_size = dp
      endif
!
      the_file%name_ = name_
!
      the_file%access_ = 'direct'
      the_file%format_ = 'unformatted'
!
      if (rec_dim .le. 0) then
         call output%error_msg("Record dimension less than zero for file "//name_)
      endif
!
      the_file%record_dim = rec_dim
      the_file%record_length = rec_dim*the_file%word_size
!
      the_file%is_open = .false.
      the_file%unit = -1
!
   end function new_direct_file
!
!
   subroutine open_direct_file(the_file, file_action)
!!
!!    Open eT direct file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file)                     :: the_file
      character(len=*), optional, intent(in) :: file_action
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if(present(file_action)) then
         the_file%action_ = trim(file_action)
      else
         the_file%action_ = 'readwrite'
      endif 
!
      if (the_file%is_open) then
!
         call output%error_msg(trim(the_file%name_)//' is already open')
!
      endif
!
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
           action=the_file%action_, recl=the_file%record_length, status='unknown', &
           form=the_file%format_, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not open eT direct file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%is_open = .true.
!
   end subroutine open_direct_file
!
!
   subroutine write_r_direct_file(the_file, scalar, record)
!!
!!    Direct file write, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      real(dp), intent(in) :: scalar
      integer, intent(in)  :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_r_direct_file
!
!
   subroutine write_r_1_direct_file(the_file, array, record)
!!
!!    Direct file write, real(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      real(dp), dimension(the_file%record_dim), intent(in) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_r_1_direct_file
!
!
   subroutine write_r_2_direct_file(the_file, array, record)
!!
!!    Direct file write, real(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      real(dp), dimension(:,:), intent(in)   :: array
      integer, intent(in)                    :: record
!
      call the_file%write_r_1_direct_file(array, record)
!
   end subroutine write_r_2_direct_file
!
!
   subroutine write_r_3_direct_file(the_file, array, record)
!!
!!    Direct file write, real(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      real(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                    :: record
!
      call the_file%write_r_1_direct_file(array, record)
!
   end subroutine write_r_3_direct_file
!
!
   subroutine write_r_4_direct_file(the_file, array, record)
!!
!!    Direct file write, real(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:,:), intent(in)  :: array
      integer, intent(in)                       :: record
!
      call the_file%write_r_1_direct_file(array, record)
!
   end subroutine write_r_4_direct_file
!
!
   subroutine write_c_direct_file(the_file, scalar, record)
!!
!!    Direct file write, complex(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      complex(dp), intent(in) :: scalar
      integer, intent(in)  :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_c_direct_file
!
!
   subroutine write_c_1_direct_file(the_file, array, record)
!!
!!    Direct file write, complex(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      complex(dp), dimension(the_file%record_dim), intent(in) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_c_1_direct_file
!
!
   subroutine write_c_2_direct_file(the_file, array, record)
!!
!!    Direct file write, complex(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)          :: the_file
      complex(dp), dimension(:,:), intent(in) :: array
      integer, intent(in)                     :: record
!
      call the_file%write_c_1_direct_file(array, record)
!
   end subroutine write_c_2_direct_file
!
!
   subroutine write_c_3_direct_file(the_file, array, record)
!!
!!    Direct file write, complex(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      complex(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                       :: record
!
      call the_file%write_c_1_direct_file(array, record)
!
   end subroutine write_c_3_direct_file
!
!
   subroutine write_c_4_direct_file(the_file, array, record)
!!
!!    Direct file write, complex(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)              :: the_file
      complex(dp), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in)                         :: record
!
      call the_file%write_c_1_direct_file(array, record)
!
   end subroutine write_c_4_direct_file
!
!
   subroutine write_i_direct_file(the_file, scalar, record)
!!
!!    Direct file write, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      integer, intent(in) :: scalar
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_direct_file
!
!
   subroutine write_i_1_direct_file(the_file, array, record)
!!
!!    Direct file write, integer 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      integer, dimension(the_file%record_dim), intent(in) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_1_direct_file
!
!
   subroutine write_i_2_direct_file(the_file, array, record)
!!
!!    Direct file write, integer 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)      :: the_file
      integer, dimension(:,:), intent(in) :: array
      integer, intent(in)                 :: record
!
      call the_file%write_i_1_direct_file(array, record)
!
   end subroutine write_i_2_direct_file
!
!
   subroutine write_i_3_direct_file(the_file, array, record)
!!
!!    Direct file write, integer 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      integer, dimension(:,:,:), intent(in)  :: array
      integer, intent(in)                    :: record
!
      call the_file%write_i_1_direct_file(array, record)
!
   end subroutine write_i_3_direct_file
!
!
   subroutine write_i_4_direct_file(the_file, array, record)
!!
!!    Direct file write, integer 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      integer, dimension(:,:,:,:), intent(in)   :: array
      integer, intent(in)                       :: record
!
      call the_file%write_i_1_direct_file(array, record)
!
   end subroutine write_i_4_direct_file
!
!
   subroutine read_r_direct_file(the_file, scalar, record)
!!
!!    Direct file read, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      real(dp), intent(out)   :: scalar
      integer, intent(in)     :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_r_direct_file
!
!
   subroutine read_r_1_direct_file(the_file, array, record)
!!
!!    Direct file read, real(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      real(dp), dimension(the_file%record_dim), intent(out) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_r_1_direct_file
!
!
   subroutine read_r_2_direct_file(the_file, array, record)
!!
!!    Direct file read, real(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      real(dp), dimension(:,:), intent(out)  :: array
      integer, intent(in)                    :: record
!
      call the_file%read_r_1_direct_file(array,record)
!
   end subroutine read_r_2_direct_file
!
!
   subroutine read_r_3_direct_file(the_file, array, record)
!!
!!    Direct file read, real(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:), intent(out)   :: array
      integer, intent(in)                       :: record
!
      call the_file%read_r_1_direct_file(array,record)
!
   end subroutine read_r_3_direct_file
!
!
   subroutine read_r_4_direct_file(the_file, array, record)
!!
!!    Direct file read, real(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                       :: record
!
      call the_file%read_r_1_direct_file(array,record)
!
   end subroutine read_r_4_direct_file
!
!
   subroutine read_c_direct_file(the_file, scalar, record)
!!
!!    Direct file read, complex(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      complex(dp), intent(out) :: scalar
      integer, intent(in)      :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_c_direct_file
!
!
   subroutine read_c_1_direct_file(the_file, array, record)
!!
!!    Direct file read, complex(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      complex(dp), dimension(the_file%record_dim), intent(out) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_c_1_direct_file
!
!
   subroutine read_c_2_direct_file(the_file, array, record)
!!
!!    Direct file read, complex(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)           :: the_file
      complex(dp), dimension(:,:), intent(out) :: array
      integer, intent(in)                      :: record
!
      call the_file%read_c_1_direct_file(array,record)
!
   end subroutine read_c_2_direct_file
!
!
   subroutine read_c_3_direct_file(the_file, array, record)
!!
!!    Direct file read, complex(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)             :: the_file
      complex(dp), dimension(:,:,:), intent(out) :: array
      integer, intent(in)                        :: record
!
      call the_file%read_c_1_direct_file(array,record)
!
   end subroutine read_c_3_direct_file
!
!
   subroutine read_c_4_direct_file(the_file, array, record)
!!
!!    Direct file read, complex(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)               :: the_file
      complex(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                          :: record
!
      call the_file%read_c_1_direct_file(array,record)
!
   end subroutine read_c_4_direct_file
!
!
   subroutine read_i_direct_file(the_file, scalar, record)
!!
!!    Direct file read, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      integer, intent(out) :: scalar
      integer, intent(in)  :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_direct_file
!
!
   subroutine read_i_1_direct_file(the_file, array, record)
!!
!!    Direct file read, integer 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      integer, dimension(the_file%record_dim), intent(out) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_1_direct_file
!
!
   subroutine read_i_2_direct_file(the_file, array, record)
!!
!!    Direct file read, integer 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      integer, dimension(:,:), intent(out)   :: array
      integer, intent(in)                    :: record
!
      call the_file%read_i_1_direct_file(array, record)
!
   end subroutine read_i_2_direct_file
!
!
   subroutine read_i_3_direct_file(the_file, array, record)
!!
!!    Direct file read, integer 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      integer, dimension(:,:,:), intent(out) :: array
      integer, intent(in)                    :: record
!
      call the_file%read_i_1_direct_file(array, record)
!
   end subroutine read_i_3_direct_file
!
!
   subroutine read_i_4_direct_file(the_file, array, record)
!!
!!    Direct file read, integer 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      integer, dimension(:,:,:,:), intent(out)  :: array
      integer, intent(in)                       :: record
!
      call the_file%read_i_1_direct_file(array, record)
!
   end subroutine read_i_4_direct_file
!
!
   subroutine write_chimera_direct_file(the_file, n_z, n_y, n_x, &
                                      & min_z, max_z, min_y, max_y, min_x, max_x, &
                                      & vector)
!!
!!    Write chimera
!!    Written by Rolf H. Myhre, September 2019
!!
!!    Specialized routine to write density grids 
!!    in the plt format readable by chimera
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
!
      integer, intent(in)  :: n_z, n_y, n_x !Number of grid points in each direction
      real(dp), intent(in) :: min_z, max_z, min_y, max_y, min_x, max_x !Min and max in each direction
!
      real(kind=4), dimension(:), intent(in) :: vector
!
      integer(i6) :: int1, int2
      integer(i6) :: n_z_sp, n_y_sp, n_x_sp !Single precision 
!
      real(kind=4) :: min_z_sp, max_z_sp, min_y_sp, max_y_sp, min_x_sp, max_x_sp !Single precision
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      int1 = 3 !Required integer that must always be 3
      int2 = 200 !Required integer that can be anything
!
      n_z_sp = int(n_z,i6)
      n_y_sp = int(n_y,i6)
      n_x_sp = int(n_x,i6)
!
      min_z_sp = real(min_z,4)
      max_z_sp = real(max_z,4)
      min_y_sp = real(min_y,4)
      max_y_sp = real(max_y,4)
      min_x_sp = real(min_x,4)
      max_x_sp = real(max_x,4)
!
      write(the_file%unit, rec=1, &
            iostat=io_error, iomsg=io_msg) int1, &
                                           int2, &
                                           n_z_sp, n_y_sp, n_x_sp, &
                                           min_z_sp, max_z_sp, &
                                           min_y_sp, max_y_sp, &
                                           min_x_sp, max_x_sp, &
                                           vector
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_chimera_direct_file
!
!
   subroutine read_chimera_direct_file(the_file, n_z, n_y, n_x, &
                                      & min_z, max_z, min_y, max_y, min_x, max_x, &
                                      & vector)
!!
!!    Read chimera
!!    Written by Andreas Skeidsvoll, Oct 2019
!!
!!    Specialized routine to read density grids 
!!    in the plt format readable by chimera
!!    
!!    Adapted from write_chimera_direct_file
!!    Written by Rolf H. Myhre, September 2019
!!
!!    Changed to from write to read
!!
      implicit none
!
      class(direct_file), intent(inout) :: the_file
!
      integer, intent(out)  :: n_z, n_y, n_x !Number of grid points in each direction
      real(dp), intent(out) :: min_z, max_z, min_y, max_y, min_x, max_x !Min and max in each direction
!
      real(kind=4), dimension(:), intent(out) :: vector
!
      integer(i6) :: int1, int2
      integer(i6) :: n_z_sp, n_y_sp, n_x_sp !Single precision 
!
      real(kind=4) :: min_z_sp, max_z_sp, min_y_sp, max_y_sp, min_x_sp, max_x_sp !Single precision
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      int1 = 3 !Required integer that must always be 3
      int2 = 200 !Required integer that can be anything
!
      read(the_file%unit, rec=1, &
            iostat=io_error, iomsg=io_msg) int1, &
                                           int2, &
                                           n_z_sp, n_y_sp, n_x_sp, &
                                           min_z_sp, max_z_sp, &
                                           min_y_sp, max_y_sp, &
                                           min_x_sp, max_x_sp, &
                                           vector
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      n_z = int(n_z_sp)
      n_y = int(n_y_sp)
      n_x = int(n_x_sp)
!
      min_z = real(min_z_sp)
      max_z = real(max_z_sp)
      min_y = real(min_y_sp)
      max_y = real(max_y_sp)
      min_x = real(min_x_sp)
      max_x = real(max_x_sp)
!
   end subroutine read_chimera_direct_file
!
!
   subroutine destructor(the_file)
!!
!!    Destructor 
!!    Written by Rolf H. Myhre, Sep 2019 
!!
      implicit none 
!
      type(direct_file) :: the_file
!
      if (the_file%is_open) then
         print *, 'Error in '//the_file%name_
         print *, 'The file is open and out of scope'
         stop
      endif
!
   end subroutine destructor
!
!
end module direct_file_class
