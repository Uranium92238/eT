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
module citation_printer_class
!
!!
!!    Citation printer
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Class for a citation printer.
!!
!!    Usage:
!!
!!       citation = new_citation(...)
!!       call citation_printer%add(citation)
!!
!!       Add more citations if you wish...
!!
!!       call citation_printer%print_(file)
!!
!
   use eT_references_class, only: eT_references
   use output_file_class,   only: output_file
!
   implicit none
!
   type :: citation_printer
!
      integer, private :: n_citations
!
      character(len=40), dimension(50), private :: implementation_labels
!
      type(eT_references), private :: references
!
      logical, private :: print_apa ! If true, prints full APA-style references;
                                    ! Otherwise, only the DOI is printed.
!
   contains
!
      procedure, public :: add &
                        => add_citation_printer
!
      procedure, public :: print_ &
                        => print_citation_printer
!
   end type citation_printer
!
   interface citation_printer
!
      procedure :: new_citation_printer
!
   end interface citation_printer
!
!  Main citation object
!
   type(citation_printer) :: eT_citations
!
contains
!
!
   function new_citation_printer(input) result(this)
!!
!!    New citation printer
!!    Written by Eirik F. Kjønstad, 2021
!!
      use input_file_class, only: input_file
!
      implicit none
!
      type(citation_printer) :: this
!
      class(input_file), intent(in) :: input 
!
      this%references = eT_references()
!
      this%implementation_labels = ''
!
!     Initially, there is just the eT reference
      this%n_citations = 1
      this%implementation_labels(1) = 'eT'
!
      this%print_apa = input%is_keyword_present('full references', 'print')
!
   end function new_citation_printer
!
!
   subroutine add_citation_printer(this, label)
!!
!!    Add
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Adds a citation.
!!
      implicit none
!
      class(citation_printer), intent(inout) :: this
!
      character(len=*), intent(in) :: label
!
      integer :: l
      logical :: new_label
!
      new_label = .true.
!
      do l = 1, this%n_citations
!
         if (trim(label) == trim(this%implementation_labels(l))) then
            new_label = .false.
            exit
         end if
!
      end do
!
      if (new_label) then
         this%n_citations = this%n_citations + 1
         this%implementation_labels(this%n_citations) = label
      end if
!
   end subroutine add_citation_printer
!
!
   subroutine print_citation_printer(this, file)
!!
!!    Print
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Prints all citations to the output file named 'file'
!!
      implicit none
!
      class(citation_printer), intent(inout) :: this
!
      class(output_file), intent(inout) :: file
!
      integer :: k
!
      character(len=:), allocatable :: doi_citation
      character(len=600) :: apa_citation
!
      call file%printf('m', '- Implementation references:', fs='(/t3,a/)')
!
      do k = 1, this%n_citations
!
         doi_citation = this%references%get_doi_citation(this%implementation_labels(k))
!
         call file%printf('m', doi_citation, fs='(t6,a)')
!
         if (this%print_apa) then
!
            apa_citation = this%references%get_apa_citation(this%implementation_labels(k))
!
            call file%printf('m', apa_citation, ffs='(/t9,a)', fs='(t9,a)')
!
            if (k .lt. this%n_citations) call file%newline('m')
!
         endif
!
      enddo
!
   end subroutine print_citation_printer
!
!
end module citation_printer_class
