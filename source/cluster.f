c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

c     contains routines :
c     cofmass - sets origin of coordinate system to centre of mass of cluster
c     rotate - rotates cluster about cofmass by randomly generated plane

C Subroutine cofmass
C Moves cluster to centre of mass

      subroutine cofmass(natoms,off)
      implicit none

C Inherited integer values
      integer natoms
C Internal integer values
      integer i
C Inherited real values
      real*8 off(3*natoms)
C Internal real values
      real*8 cofmx
      real*8 cofmy
      real*8 cofmz
      real*8 rnat
      
      cofmx = 0.d0
      cofmy = 0.d0
      cofmz = 0.d0
      
      do i=1,natoms
C Sum all coordinate totals
         cofmx = cofmx + off(i)
C      end do
C      do i=natoms+1,2*natoms
         cofmy = cofmy + off(natoms+i)
C      end do
C      do i=2*natoms+1,3*natoms
         cofmz = cofmz + off((2*natoms)+i)
      end do

      rnat = 1/dble(natoms)	! Reciprocal of number of atoms

      cofmx = cofmx*rnat	! Scale x coordinate
      cofmy = cofmy*rnat 	! Scale y coordinate
      cofmz = cofmz*rnat	! Scale z coordinate

      do i=1,natoms
         off(i) = off(i) - cofmx
C      end do
C      do i=natoms+1,2*natoms
         off(natoms+i) = off(natoms+i) - cofmy
C      end do
C      do i=2*natoms+1,3*natoms
         off((2*natoms)+i) = off((2*natoms)+i) - cofmz
      end do

      return
      end

C Subroutine rotate
c     Rotates the two clusters in the array pair by
c     different randomly generated values of theta and phi

      SUBROUTINE rotate(natoms,off,idum,pi)
      implicit none

C Inherited integer values
      integer natoms
      integer idum	! Random number seed
C Internal integer values
      integer i
C Inherited real values
      real*8 off(3*natoms)
      real*8 pi
C Internal real values
      real*8 rot11
      real*8 rot12
      real*8 rot13
      real*8 rot21
      real*8 rot22
      real*8 rot23
      real*8 rot31
      real*8 rot32
      real*8 rot33
      real*8 theta
      real*8 phi
      real*8 ratomsx
      real*8 ratomsy
      real*8 ratomsz
      real*8 ran3
C Import random number generator
      external ran3

c     generate random plane

      theta = ran3(idum)*(pi/2.d0)	! Random number in range 2pi
      phi = ran3(idum)*pi		! Random number in range pi

c     Sets up rotation matrix
C As standard for any x,y,z rotation matrix

      rot11 = dble(cos(phi))
      rot12 = 0.d0
      rot13 = dble(-sin(phi))
      rot21 = dble(sin(theta)*sin(phi))
      rot22 = dble(cos(theta))
      rot23 = dble(sin(theta)*cos(phi))
      rot31 = dble(cos(theta)*sin(phi))
      rot32 = dble(-sin(theta))
      rot33 = dble(cos(theta)*cos(phi))

      ratomsx = 0.d0		! Set for temp variables x
      ratomsy = 0.d0		! Set for temp variables y
      ratomsz = 0.d0		! Set for temp variables z

c     rotation - no innner loops
c     first cluster in 'pair'

      do i=1,natoms
C For all atoms
         ratomsx = rot11*off(i) + rot12*off(i+natoms)
     &        + rot13*off(i+2*natoms)

         ratomsy = rot21*off(i) + rot22*off(i+natoms)
     &        + rot23*off(i+2*natoms)

         ratomsz = rot31*off(i) + rot32*off(i+natoms)
     &        + rot33*off(i+2*natoms)

C Save rotations in offspring positions
         off(i)          = ratomsx
         off(i+natoms)   = ratomsy
         off(i+2*natoms) = ratomsz
      
      end do

      return
      end

