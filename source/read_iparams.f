c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Subroutine to read in Rigid- Ion potentials
      subroutine read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &            B_bb,rho_bb,bscale,namea,nameb)
      implicit none

C     Integer Value for file length
      integer ifilelen
 
C     Real values for variables
      real*8 r_ab
      real*8 qa
      real*8 qb
      real*8 B_ab
      real*8 rho_ab
      real*8 B_bb
      real*8 rho_bb
      real*8 bscale

C     Strings for variables
      character*(*)inputfile
      character*(*)namea
      character*(*)nameb

C     Open input file
      open(unit=2,status='old',file=inputfile(1:ifilelen))

C     Read in variables in order
      read(2,*)r_ab		! Bond distance
C     Atom types on same line
      read(2,'(a,a)')namea,nameb	! Atom types
      read(2,*)qa,qb		!
      read(2,*)B_ab,rho_ab	!
      read(2,*)B_bb,rho_bb	!

C     Scaling factor
      read(2,*)bscale

C     Close file and exit
      close(2)

      return
      end
