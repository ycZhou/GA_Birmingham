c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Subroutine to read in cutoff variables
      subroutine read_cutoff(cutoff_start,cutoff_end,a3a,a4a,a5a,
     & x3x,x4x,x5x)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C     Local integer values needed
      integer nspec	! Number of atoms types
      integer i		! Current array position count

C     Arrays of real numbers, global variables
      real*8 cutoff_start(nspecmax2)		! Cutoff start distance
      real*8 cutoff_end(nspecmax2)		! Cutoff end distance
      real*8 a3a(nspecmax2)			! a3 Coefficient
      real*8 a4a(nspecmax2)			! a4 Coefficient
      real*8 a5a(nspecmax2)			! a5 Coefficient
      real*8 x3x(nspecmax2)			! x3 Coefficient
      real*8 x4x(nspecmax2)			! x4 Coefficient
      real*8 x5x(nspecmax2)			! x5 Coefficient

c     Read in parameters
C     Firstly we'll open up the cutoff_parameters file.      
      open(unit=2,file='cutoff_parameters',status='old')

C     Read in number of elements (Currently configured up to bimetallic
      read(2,*)nspec

C     Set initial array positioning
      i=1

C     Read in opening details
C     There are no quality checks on this but we must assume the user will
C     get this information correct.
      read(2,*)cutoff_start(i)
      read(2,*)cutoff_end(i)
      read(2,*)a3a(i)
      read(2,*)a4a(i)
      read(2,*)a5a(i)
      read(2,*)x3x(i)
      read(2,*)x4x(i)
      read(2,*)x5x(i)

C     If number of elements is 2, then add hetrogeneous bonds and b-b type bonds to array
C     In the order a-b, b-a, bb
      if (nspec .eq. 2) then

C         Firstly a-b
          i=i+1
          read(2,*)cutoff_start(i)
          read(2,*)cutoff_end(i)
          read(2,*)a3a(i)
          read(2,*)a4a(i)
          read(2,*)a5a(i)
          read(2,*)x3x(i)
          read(2,*)x4x(i)
          read(2,*)x5x(i)

C         Then duplicate as b-a. This are not read from file
          i=i+1
          cutoff_start(i)=cutoff_start(i-1)
          cutoff_end(i)=cutoff_end(i-1)
          a3a(i)=a3a(i-1)
          a4a(i)=a4a(i-1)
          a5a(i)=a5a(i-1)
          x3x(i)=x3x(i-1)
          x4x(i)=x4x(i-1)
          x5x(i)=x5x(i-1)
                                                       
C         Then read in b-b
          i=i+1
          read(2,*)cutoff_start(i)
          read(2,*)cutoff_end(i)
          read(2,*)a3a(i)
          read(2,*)a4a(i)
          read(2,*)a5a(i)
          read(2,*)x3x(i)
          read(2,*)x4x(i)
          read(2,*)x5x(i)

      end if

C     Close file and exit
      close(unit=2)

      return
      end
