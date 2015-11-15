c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Subroutine to read morse potential
      subroutine read_morse(inputfile,ifilelen,alpha,bscale,el)
      implicit none

C     Integer for local variable
      integer ifilelen	! Filename length

C     Morse parameters
      real*8 alpha
      real*8 bscale
      character*(*)el
      character*(*)inputfile

c     ... read in potential parameters
C     Open file
      open(2,status='old',file=inputfile(1:ifilelen))

C     Read variables
      read(2,*)alpha
      read(2,'(a,f5.4)')el,bscale	! Read in element type and scalor
      
C     Close file and return
      close(2)

      return
      end
