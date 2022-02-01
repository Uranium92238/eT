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
module accelerator_factory_class
!
!!
!! Accelerator factory class
!! Written by Sarai D. Folkestad, 2020
!!
!! Constructs an accelerator for solving equations
!!
!!    Possible accelerators:
!!    - DIIS (or CROP)
!!    - None
!!
!
   use parameters
!
   use accelerator_tool_class,      only: accelerator_tool
   use null_accelerator_tool_class, only: null_accelerator_tool
   use diis_accelerator_tool_class, only: diis_accelerator_tool
   use diis_tool_class,             only: diis_tool
!
   implicit none
!
   type :: accelerator_factory
!
      character(len=200), private :: input_section
!
   contains
!
      procedure :: create => create_accelerator_factory
!
      procedure, private :: read_diis_settings &
                         => read_diis_settings_accelerator_factory
!
   end type accelerator_factory
!
   interface accelerator_factory
!
      procedure :: new_accelerator_factory
!
   end interface accelerator_factory
!
contains
!
!
   function new_accelerator_factory(input_section) result(this)
!!
!!    New accelerator factory
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      type(accelerator_factory)  :: this
      character(len=*)           :: input_section
!
      this%input_section = input_section
!
   end function new_accelerator_factory
!
   function create_accelerator_factory(this,&
                                       accelerator_type, &
                                       x_dimension,      &
                                       e_dimension) result(accelerator)
!!
!!    Create
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Constructs the accelerator
!!
      use global_out, only: output
!
      implicit none
!
      class(accelerator_factory), intent(in) :: this
      class(accelerator_tool), allocatable   :: accelerator
!
      character(len=*), intent(in)           :: accelerator_type
      integer, intent(in)                    :: x_dimension
      integer, intent(in)                    :: e_dimension
!
      logical            :: crop
      logical            :: records_in_memory
      integer            :: diis_dimension
      character(len=200) :: diis_name
!
      if (trim(accelerator_type) == 'none') then
!
         accelerator = null_accelerator_tool(x_dimension, e_dimension)
!
      elseif (trim(accelerator_type) == 'diis') then
!
         diis_name         = trim(this%input_section)
         diis_dimension    = 8
         crop              = .false.
         records_in_memory = .true.
!
         call this%read_diis_settings(diis_dimension, crop, records_in_memory)
!
         accelerator = diis_accelerator_tool(x_dimension, e_dimension,                 &
                                             diis_tool(name_        = diis_name,       &
                                                       n_parameters = x_dimension,     &
                                                       n_equations  = e_dimension,     &
                                                       dimension_   = diis_dimension,  &
                                                       crop         = crop,            &
                                                       records_in_memory = records_in_memory))
!
      else
!
         call output%error_msg('in accelerator factory. Did not recognize accelerator type')
!
      endif
!
   end function create_accelerator_factory
!
!
   subroutine read_diis_settings_accelerator_factory(this, diis_dimension, crop, records_in_memory)
!!
!!    Read DIIS settings
!!    Written by Sarai D. Folkestad
!!
      use global_in, only: input
!
      implicit none
!
      class(accelerator_factory),   intent(in)    :: this
      integer,                      intent(inout) :: diis_dimension
      logical,                      intent(inout) :: crop
      logical,                      intent(inout) :: records_in_memory
!
      character(len=200) :: storage
!
      call input%get_keyword('diis dimension', trim(this%input_section), diis_dimension)
!
      crop = input%is_keyword_present('crop', trim(this%input_section))
!
      call input%get_keyword('storage', trim(this%input_section), storage)
!
      if (trim(storage) == 'memory') then
!
         records_in_memory = .true.
!
      elseif (trim(storage) == 'disk') then
!
         records_in_memory = .false.
!
      endif
!
   end subroutine read_diis_settings_accelerator_factory
!
!
end module accelerator_factory_class
