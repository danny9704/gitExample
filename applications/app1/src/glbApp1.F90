module convert_library
   implicit none
   
contains
   subroutine convert_pat_to_hst()
      implicit none
      integer :: iounit_in, iounit_out, i
      character(len=256) :: line
      character(len=256), parameter :: input_file = 'Semi.pat'
      character(len=256), parameter :: output_file = 'Semi.hst'

      open(unit=iounit_in, file=input_file, status='old', action='read')
      open(unit=iounit_out, file=output_file, status='replace', action='write')

      do
         read(iounit_in, '(A)', iostat=i) line
         if (i /= 0) exit
         write(iounit_out, '(A)') trim(line)
      end do

      close(iounit_in)
      close(iounit_out)
   end subroutine convert_pat_to_hst
end module convert_library
