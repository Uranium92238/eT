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
module citation_class
!
!!
!!    Citation
!!    Written by Eirik F. Kjønstad, 2021
!!
!
   implicit none
!
   type :: citation
!
      character(len=:), allocatable, private :: journal
      character(len=:), allocatable, private :: title_
      character(len=:), allocatable, private :: volume
      character(len=:), allocatable, private :: issue
      character(len=:), allocatable, private :: pages
      character(len=:), allocatable, private :: year
      character(len=:), allocatable, private :: doi
      character(len=:), allocatable, private :: implementation
!
!     With GCC-9 we could use len=: and
!     allocate(this%authors, source=authors) in the constructor
      character(len=25), dimension(:), allocatable, private :: authors
!
   contains
!
      procedure, public :: get_implementation
      procedure, public :: get_doi_citation_string
      procedure, public :: get_apa_citation_string
!
   end type citation
!
   interface citation
!
      procedure :: new_citation
!
   end interface citation
!
contains
!
!
   function new_citation(implementation, journal, title_, volume, issue, &
                         pages, year, doi, authors) result(this)
!!
!!    New citation printer
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Implementation: describes what method/algorithm/implementation the
!!                    citation refers to.
!!
      implicit none
!
      type(citation) :: this
!
      character(len=*), intent(in) :: implementation
      character(len=*), intent(in) :: journal
      character(len=*), intent(in) :: title_
      character(len=*), intent(in) :: volume
      character(len=*), intent(in) :: issue
      character(len=*), intent(in) :: pages
      character(len=*), intent(in) :: year
      character(len=*), intent(in) :: doi
!
      character(len=*), dimension(:), intent(in) :: authors
!
      this%implementation = implementation
      this%journal        = journal
      this%title_         = title_
      this%volume         = volume
      this%issue          = issue
      this%pages          = pages
      this%year           = year
      this%doi            = doi
!
      allocate(this%authors(size(authors)))
      this%authors = authors
!
   end function new_citation
!
!
   function get_implementation(this) result(implementation)
!!
!!    Get implementation
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(citation), intent(in) :: this
!
      character(len=:), allocatable :: implementation
!
      implementation = trim(this%implementation)
!
   end function get_implementation
!
!
   function get_doi_citation_string(this) result(doi_citation)
!!
!!    Get DOI citation string
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(citation), intent(in) :: this
!
      character(len=:), allocatable :: doi_citation
!
      doi_citation = this%implementation // ": https://doi.org/" // this%doi
!
   end function get_doi_citation_string
!
!
   function get_apa_citation_string(this) result(string)
!!
!!    Get APA citation string
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Makes a citation string from the reference information: author, year, journal, etc.
!!
!!    The string resembles APA-style, e.g.:
!!
!!       Sarai D. Folkestad, Eirik F. Kjønstad, and Henrik Koch (2019). An
!!       efficient algorithm for Cholesky decomposition of electron repulsion
!!       integrals. J. Chem. Phys. 150(19), 194112.
!!
      use global_out, only: output
!
      implicit none
!
      class(citation), intent(in) :: this
!
      character(len=600) :: string
!
      integer :: author, n_authors
!
      n_authors = size(this%authors)
!
      if (n_authors             .eq. 0 .or.  &
          len_trim(this%title_) .eq. 0 .or.  &
          len_trim(this%year)   .eq. 0) then
!
         call output%error_msg('a citation must have both authors, title, and year.')
!
      endif
!
!     Authors
!
      if (n_authors == 1) then
!
         string = trim(this%authors(1))
!
      else
!
         string = ''
!
         do author = 1, n_authors - 1
!
            string = trim(string) // ' ' // trim(this%authors(author)) // ','
!
         enddo
!
         write(string, '(3a)') trim(string), ' and ', trim(this%authors(n_authors))
!
      end if
!
!     Year
!
      write(string, '(4a)') trim(string), ' (', trim(this%year), ')'
!
!     Title, journal, volume, issue, pages
!
      write(string, '(3a)') trim(string), '. ', trim(this%title_)
!
      if (len_trim(this%journal) .gt. 0) then
!
            write(string, '(3a)') trim(string), '. ', trim(this%journal)
!
      endif
!
      if (len_trim(this%volume) .gt. 0) then
!
            if (len_trim(this%journal) .eq. 0) &
               call output%error_msg('a citation with a volume must also have a journal.')
!
            write(string, '(3a)') trim(string), ' ', trim(this%volume)
!
      endif
!
      if (len_trim(this%issue) .gt. 0) then
!
            if (len_trim(this%journal) .eq. 0) &
               call output%error_msg('a citation with an issue must also have a journal.')
!
            write(string, '(4a)') trim(string), ' (', trim(this%issue), ')'
!
      endif
!
      if (len_trim(this%pages) .gt. 0) then
!
            if (len_trim(this%journal) .eq. 0) &
               call output%error_msg('a citation with pages must also have a journal.')
!
            write(string, '(3a)') trim(string), ', ', trim(this%pages)
!
      endif
!
      string = trim(string) // '.'
!
      string = adjustl(string)
!
   end function get_apa_citation_string
!
!
end module citation_class
