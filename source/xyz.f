      subroutine xyz(natoms,pfunc,off,qoff,bscale,suff,int,flen,
     &     namea,nameb,el,fname,nels)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

      integer natoms,pfunc,int,flen,nels
      real*8 off(3*natoms+1),qoff(natoms),bscale

      character*(*)namea,nameb,el
      character*(fnamelen)fname

      integer i,j,len,in
      character*20 suff,str

c     writes out coordinates to file 'fname'.xyz in xyz format

c     copy suffix

      if (suff(1:1) .eq. 'i') then
         in = 1
      else
         in = int
      end if
      
c     determines length of suffix and adds integer to end

      call string(in,suff,str,len)

      if (str(1:len) .eq. 'none00' ) then
         open(3,status='unknown',file=fname(1:flen)//'.xyz')
      else
         open(3,status='unknown',
     &        file=fname(1:flen)//str(1:len)//'.xyz')
      end if

      write(3,*)natoms
      write(3,'(f13.6,2x,f7.4)')off(3*natoms+1),bscale
      if (nels .eq. 2) then
         if (pfunc .eq. 3 .or. pfunc .eq. 5) then
            do j=1,natoms
               if (qoff(j) .gt. 0.d0) then
                  write(3,'(a,1x,3f13.8)')namea,bscale*off(j),
     &                 bscale*off(j+natoms),bscale*off(j+2*natoms)
               else
                  write(3,'(a,1x,3f13.8)')nameb,bscale*off(j),
     &                 bscale*off(j+natoms),bscale*off(j+2*natoms)
               end if
            end do
         else if (pfunc .eq. 4 .or. pfunc .eq. 6) then
            do j=1,natoms
               if (qoff(j) .eq. 1.0d0) then
                  write(3,'(a,1x,3f13.8)')namea,bscale*off(j),
     &                 bscale*off(j+natoms),bscale*off(j+2*natoms)
               else
                  write(3,'(a,1x,3f13.8)')nameb,bscale*off(j),
     &                 bscale*off(j+natoms),bscale*off(j+2*natoms)
               end if
            end do
         end if
      else
         if (pfunc .eq. 4 .or. pfunc .eq. 6 ) then
            do j=1,natoms
               write(3,'(a,1x,3f13.8)')namea,bscale*off(j),
     &              bscale*off(j+natoms),bscale*off(j+2*natoms)
            end do
         else 
            do j=1,natoms
               write(3,'(a,1x,3f13.8)')el,bscale*off(j),
     &              bscale*off(j+natoms),bscale*off(j+2*natoms) 
            end do
         end if
      end if
      
      close(3)

      return
      end
