c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Subrotuine to read in gupta potential
      subroutine read_gupta(inputfile,ifilelen,nspec,aij,vij,pij,qij,r0,
     &     bscale,namea,nameb)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C     Local integer values
      integer ifilelen	! File name length (?neccessary?)
      integer nspec	! Number of elements
      integer i		! Current array location reference

C     Arrays of potential variables
      real*8 aij(nspecmax2)	! A
      real*8 vij(nspecmax2)	! zi
      real*8 pij(nspecmax2)	! p
      real*8 qij(nspecmax2)	! q
      real*8 r0(nspecmax2)	! Bond lengths
      real*8 bscale		!

C     String variables
      character*(*)inputfile	! Input filename
      character*(*)namea	! Element a
      character*(*)nameb	! Element b (if needed)

C      (Original comments)
C      integer i,j,ii
C      Looks like these come from trimetallic setup

C     read in potential
C     Open file     
      open(unit=2,file=inputfile(1:ifilelen),status='old')

C     Read number of elements
      read(2,*)nspec

C     Set current location in array for a-a interactions
      i=1

C     Read in variables
      read(2,*)aij(i)	! A
      read(2,*)vij(i)	! zi
      read(2,*)pij(i)	! p
      read(2,*)qij(i)	! q
      read(2,*)r0(i)	! r0
 
C     If number of elements is two...
      if (nspec .eq. 2) then
          
C         Increase array pointer
          i=i+1
C         And read in next variables for a-b interactions
          read(2,*)aij(i)
          read(2,*)vij(i)
          read(2,*)pij(i)
          read(2,*)qij(i)
          read(2,*)r0(i)

C         Duplicate a-b for b-a interactions
          i=i+1
          aij(i)=aij(i-1)
          vij(i)=vij(i-1)
          pij(i)=pij(i-1)
          qij(i)=qij(i-1)
          r0(i)=r0(i-1)

C         Then read in b-b interactions
          i=i+1
          read(2,*)aij(i)
          read(2,*)vij(i)
          read(2,*)pij(i)
          read(2,*)qij(i)
          read(2,*)r0(i)

      end if

C      OUTDATED CODE
C      This is from when the system was configured for trimetallics. Could be reimplemented
C      Input file would need to list interactions in forms:
C      a-a, a-b, a-c, b-a, b-b, b-c, c-a, c-b, c-c
C      Obviously some of these are duplicates, and this raises inconsistency errors
C      So would be nice to change the code to deal with duplicate variable types. 
C   
C      ii=0
C      do i=1,nspec
C         do j=1,nspec
C            ii=ii+1
C            read(2,*)aij(ii)
C            read(2,*)vij(ii)
C            read(2,*)pij(ii)
C            read(2,*)qij(ii)
C            read(2,*)r0(ii)
C         enddo
C      enddo
C

C     Read in element a
      read(2,'(a)')namea
C     If bimetallic read in element b
      if (nspec .eq. 2) then
          read(2,'(a)')nameb
      end if

C     Read in scale factor (normally set to 1)
      read(2,*)bscale

C     Close file and return
      close(unit=2)

      return
      end
