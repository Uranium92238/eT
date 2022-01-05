!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module environment_factory_class
!
!!
!!    Environment_factory class module
!!    Written by Sarai D. Folkestad, Oct 2020
!!
!
   use environment_class,                 only: environment
   use pcm_environment_class,             only: pcm_environment
   use fq_environment_class,              only: fq_environment
   use electrostatic_environment_class,   only: electrostatic_environment
!
   use global_out,                        only: output
   use global_in,                         only: input
!
   implicit none
!
   type :: environment_factory
!
      character(len=200), private :: type_ 
      character(len=200), private :: forcefield 
!
   contains
!
      procedure :: create &
                => create_environment_factory
!
      procedure, private :: read_embedding &
                         => read_embedding_environment_factory
!
   end type  environment_factory
!
   interface environment_factory
!
      procedure :: new_environment_factory
!
   end interface environment_factory
!
contains
!
!
   function new_environment_factory() result(factory)
!!
!!    New environment factory
!!    Written by Sarai D. Folkestad, Oct 2020 
!!
      implicit none
!
      type(environment_factory) :: factory
!
      call factory%read_embedding()
!
   end function new_environment_factory
!
!
   function create_environment_factory(factory) result(embedding)
!!
!!    Create
!!    Written by Sarai D. Folkestad
!!
!!    Creates an embedding object based on input information.
!!
!
      use ao_tool_class, only: ao_tool
!
      implicit none
!
      class(environment_factory), intent(in) :: factory
!
      class(environment), allocatable :: embedding
!
      if (trim(factory%type_) == 'pcm') embedding = pcm_environment()
!
      if (trim(factory%type_) == 'molecular mechanics') then
!
         if (trim(factory%forcefield) == 'non-polarizable') embedding = electrostatic_environment()
         if (trim(factory%forcefield) == 'fq')              embedding = fq_environment()
!
      endif
!
   end function create_environment_factory
!
!
   subroutine read_embedding_environment_factory(factory)
!!
!!    Read embedding 
!!    Written by Sarai D. Folkestad, Oct 2020
!!
      implicit none

      class(environment_factory), intent(inout) :: factory
!
!     Cannot perform PCM and MM at the same time
!
      if (input%is_section_present('molecular mechanics') .and. &
          input%is_section_present('pcm')) &
            call output%error_msg('PCM and QM/MM together is not currently supported.')
!
      if (.not. (input%is_section_present('molecular mechanics') .or. &
                 input%is_section_present('pcm'))) &
            call output%error_msg('embedding type must be either PCM or molecular mechanics.')
!
      if (input%is_section_present('molecular mechanics')) then
!         
         factory%type_ = 'molecular mechanics'
!
         call input%get_keyword('forcefield', 'molecular mechanics', factory%forcefield)
!
!        Check if force field is supported
!
         if (trim(factory%forcefield) .ne. 'non-polarizable' .and. &
            trim(factory%forcefield) .ne. 'fq') call output%error_msg('force field not recognized.')
!
      endif
!
      if (input%is_section_present('pcm')) factory%type_ = 'pcm'
!
   end subroutine read_embedding_environment_factory
!
!
end module environment_factory_class
