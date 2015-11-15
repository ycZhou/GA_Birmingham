      subroutine tokeniser(ttok,tlen,blen,bloc,lno,buffer)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

      character*(tokbufflen) buffer
      character*(maxtoklen) token
      character*80 line
      character*20 ifile

      integer i,j,k,last,len(maxtoklen),n,tok,loc(maxtokens),tlen,ttok,
     &     bloc(maxtokens),blen(maxtokens),lent,lno(maxtokens),ln

      real*8 chartodb,rr
      logical real,eof

      tlen = 0
      ttok = 0
      eof = .false.
      k = 0

      do i=1,maxtokens
         bloc(i) = 0
         blen(i) = 0
         lno(i)  = 0
      end do

      do while ( .not. eof )

c     read in a line (80 characters worth)

         read(2,'(80a1)',end=99) (line(i:i), i=1,80)

         do i=1,maxtoklen
            loc(i)=0
            len(i)=0
         end do

c     line number

         k = k + 1

c     find last non-blank character

         last = 80
         
         do while (ichar(line(last:last)) .le. 32
     &        .and. last .gt. 0)
            last = last -1
         end do

c     find length and location of each token
c     ?expect? tokens to be separated by spaces and/or equals sign
         
         tok = 0 
         n = 1
         
         if ( last .gt. 0 ) then
            do while (n .le. last)
               if (ichar(line(n:n)) .gt. 32 
     &              .and. ichar(line(n:n)) .ne. 61) then
                  if (ichar(line((n-1):(n-1))) .eq. 61
     &                 .or. ichar(line((n-1):(n-1))) .eq. 32) then
                     tok = tok + 1
                     loc(tok+1) = n
                  end if
                  len(tok+1) = len(tok+1) + 1
                  if (n .eq. 1) loc(tok+1) = 1
                  n = n + 1
               else if (ichar(line(n:n)) .eq. 61) then
                  tok = tok + 1
                  loc(tok+1) = n
                  len(tok+1) = 1
                  n = n + 1
               else 
                  n = n + 1
               end if
            end do

c     copy tokens into buffer
c     and set arrays with locations and lengths of tokens
            
            do i=1,tok+1
               do j=1,len(i)
             buffer((tlen+j):(tlen+j)) = line((loc(i)+j-1):(loc(i)+j-1))
               end do
c     total tokens
               ttok = ttok + 1
c     starting locations of tokens in string
               bloc(ttok) = tlen + 1
c     length of tokens
               blen(ttok) = len(i)
c     line in input file token came from
               lno(ttok) = k
c     total length of master token string
               tlen = tlen + len(i)
            end do

         end if

      end do

      close(2)

c     routine exits loop of lines in input file here
     
 99   eof = .true.

      return
      end


      integer function chartoint(cnum,len)
      implicit none
      character*(*) cnum
      integer len,no

      if (len .eq. 0) goto 99

#ifdef AIX
      read(cnum,10,err=99)no
 10   format(i<len>)
#else
      read(cnum,'(i10)',err=99)no
#endif

      chartoint = no
      return

 99   chartoint = -1
      
      return      
      end

      
      real*8 function chartodb(cnum,len)
      implicit none
      character*(*) cnum
      integer len
      real*8 db

      if (len .eq. 0) goto 99

#ifdef AIX
      read(cnum,10,err=99)db
 10   format(f<len>.0)
#else
      read(cnum,'(f13.6)',err=99)db
#endif
      
      chartodb = db
      return

 99   chartodb = -1.d0

      return
      end


      subroutine get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,
     &     n)
      implicit none

      character*(*) buffer
      character*(*) token

      integer ttok,tlen,blen(ttok),bloc(ttok),lno(ttok),lent,ln,n
 
      integer i

c     gets token n and it's length and line number

      do i=1,blen(n)
         token(i:i) = buffer((bloc(n)+i-1):(bloc(n)+i-1))
      end do
      
c     check to see if we are asked for a token beyond ttok

      lent = blen(n)
      ln = lno(n)

      return
      end
      




      
