c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Subroutine to read in MM parameters
      subroutine read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,
     &               idamp,ihard,ram,rcut2,bscale,el,title)
      implicit none

C     Integer values for reading potential
      integer ifilelen	! Filename length
      integer idamp	!
      integer ihard	!

C     Real values
      real*8 de		!
      real*8 re		!
      real*8 a2		!
      real*8 a3		!
      real*8 coeff(11)	! Coefficients Array of size 11
      real*8 ram	!
      real*8 rcut2	!
      real*8 bscale	! Scaling factor
      real*8 rcut	!

C     Characters / Strings
      character*(*)inputfile	! Inputfile name
      character*(*)el		! Element
      character*(*)title	!

c     ... read in potential parameters
C     Open file
      open(2,status='old',file=inputfile(1:ifilelen))
      
C     Read in variables
      read(2,'(a)')title
      read(2,*)de,re
      read(2,*)a2,a3

C     Read in coefficients. Set to 0 if not used
      read(2,*)coeff(1),coeff(2),coeff(3),coeff(4),coeff(5),coeff(6),
     &          coeff(7)
      read(2,*)coeff(8),coeff(9),coeff(10),coeff(11)
      read(2,*)idamp,ihard
      read(2,*)ram
      read(2,*)rcut
C     Read in element type and scale
      read(2,'(a,f5.4)')el,bscale
      
      rcut2 = rcut*rcut		! rcut squared

C     Close file and return
      close(2)

      return
      end
