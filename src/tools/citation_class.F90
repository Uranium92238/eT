!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   use global_out, only: output
   use output_file_class, only: output_file
!
   implicit none
!
   type :: citation 
!
      character(len=600), private :: journal 
      character(len=600), private :: title_ 
      character(len=600), private :: volume 
      character(len=600), private :: issue 
      character(len=600), private :: pages 
      character(len=600), private :: year 
      character(len=600), private :: doi 
      character(len=600), private :: implementation 
!
      character(len=600), dimension(:), allocatable, private :: authors  
!
   contains 
!
      procedure, public :: get_doi
      procedure, public :: get_implementation
!
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
      this%authors        = authors
!
   end function new_citation
!
!
   function get_doi(this) result(doi)
!!
!!    Get DOI
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(citation), intent(in) :: this 
!
      character(len=600) :: doi 
!
      doi = this%doi
!
   end function get_doi
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
      character(len=600) :: implementation 
!
      implementation = this%implementation
!
   end function get_implementation
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
      implicit none 
!
      class(citation), intent(in) :: this 
!
      character(len=600) :: string 
!
      character(len=600) :: and 
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
      string = ''
!
      do author = 1, n_authors - 1
!
         string = trim(string) // ' ' // trim(this%authors(author)) // ','
!
      enddo
!
      and = ''
      if (n_authors .gt. 1) and = 'and'
!
      string = trim(string) // ' ' // trim(and) // ' ' // trim(this%authors(n_authors)) 
!
!     Year
!
      string = trim(string) // ' (' // trim(this%year) // ')'
!
!     Title, journal, volume, issue, pages
!
      string = trim(string)         &
      // '. '                       &
      // trim(this%title_)          
!
      if (len_trim(this%journal) .gt. 0) then
!
            string = trim(string)         &
            // '. '                       &
            // trim(this%journal)         
!
      endif
!
      if (len_trim(this%volume) .gt. 0) then 
!
            if (len_trim(this%journal) .eq. 0) &
               call output%error_msg('a citation with a volume must also have a journal.')
!
            string = trim(string)         &
            // ' '                        &
            // trim(this%volume)          
!
      endif
!
      if (len_trim(this%issue) .gt. 0) then
!
            if (len_trim(this%journal) .eq. 0) &
               call output%error_msg('a citation with an issue must also have a journal.')
!
            string = trim(string)         &
            // ' '                        &
            // '('                        & 
            // trim(this%issue)           &
            // ')'                      
!
      endif
!
      if (len_trim(this%pages) .gt. 0) then 
!
            if (len_trim(this%journal) .eq. 0) &
               call output%error_msg('a citation with pages must also have a journal.')
!
            string = trim(string)         &
            // ', '                       & 
            // trim(this%pages)           
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
