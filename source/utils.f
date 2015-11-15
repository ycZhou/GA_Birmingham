c     useful functions
c     string - converts integer into string and appends to a character array
c     lstr - finds length of character string

      SUBROUTINE string(in,suff,str,len)
      implicit none    
      integer in,len,lstr
      character*(*) suff,str
      character*4 inch
  
      external lstr

c     find length of 'suff' and copy into 'str'

      len = lstr(suff)

      str(1:len) = suff(1:len)
   
c     copy 'in' to character buffer

      write(inch,'(i4)')in

c     format according to size of 'in'

      if (in .lt. 10) then
         str(len+1:len+1) = '0'
         str(len+2:len+2) = inch(4:4)
         len = len + 2
      else if (in .lt. 100) then
         str(len+1:len+1) = '0'
         str(len+2:len+3) = inch(3:4)
         len = len + 3
      else if (in .lt. 1000) then
         str(len+1:len+3) = inch(2:4)
         len = len + 3
      else
         str(len+1:len+4) = inch(1:4)
         len = len + 4
      end if

      return
      end


      function lstr(string)
      implicit none
      character*(*) string
      integer lstr,len

      len = 0

      do while (ichar(string((len+1):(len+1))) .gt. 32) 
         len = len+1
      end do

      lstr = len

      return
      end

      


      
