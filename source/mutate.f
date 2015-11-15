c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

c     mutating routines:
c     mutate_move - moves ~1/10 atoms to random positions
c     mutate_rotate - rotates 'half' of cluster by random angle
c     mutate_replace - generates entire new cluster
c     mutate_project - projects atoms on to surface of a sphere

C Subroutine mutate_move

      SUBROUTINE mutate_move(natoms,off,idum,range)
      implicit none

C Inherited integer values
      integer natoms
      integer idum
C Internal integer values
      integer i
      integer iatom
C Inherited real values
      real*8 off(3*natoms)
      real*8 range
C Internal real values
      real*8 range2
      real*8 ran3

c     Mutates one of the new strings made during mating
   
C Import random number generator
      external ran3

c     ** incorrect version ** atoms only put in +ve quadrant **

c     Picks a random atom from a random string and changes its coordinates

c      do i=1,natoms/10
c         iatom = ran3(idum)*natoms+1
c         off(iatom) = ran3(idum)*range
c         off(iatom+natoms) = ran3(idum)*range
c         off(iatom+2*natoms) = ran3(idum)*range
c      end do

c     ** correct version **

      range2 = range/2.d0

      do i=1,natoms/3
C For 1/3 of all atoms
         iatom = ran3(idum)*natoms+1		! Pick random atom
         off(iatom) = ran3(idum)*range2			! Set random x coordinate
         off(iatom+natoms) = ran3(idum)*range2		! Set random y coordinate 
         off(iatom+2*natoms) = ran3(idum)*range2	! Set random z coordinate
         if (ran3(idum) .gt. 0.5d0) 
     &        off(iatom) = -off(iatom)				! Make x coordinate negatie with 50% probability
         if (ran3(idum) .gt. 0.5d0) 
     &        off(iatom+natoms) = -off(iatom+natoms) 		! Make y coordinate negatie with 50% probability
         if (ran3(idum) .gt. 0.5d0) 
     &        off(iatom+2*natoms) = -off(iatom+2*natoms)	! Make z coordinate negatie with 50% probability
      end do

      return	! Return
      end

C Subroutine mutate_rotate
C Mutates a cluster my rotating 'half'

      SUBROUTINE mutate_rotate(natoms,off,qoff,idum,pi)
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer idum	! Random seed
C Internal integer values
      integer i
C Inherited Real Values
      real*8 off(3*natoms)	! Coordiantes of offspring
      real*8 qoff(natoms)	! Charge of offspring
      real*8 pi			! pi
C Internal real variables
      real*8 ran3
      real*8 x
      real*8 y
      real*8 theta
      real*8 rot11
      real*8 rot12
      real*8 rot21
      real*8 rot22
C Import random number      
      external ran3
      
C (Original Comment)
c     zsort cluster
C Sort cluster by z coordinates

      call z_heapsort(natoms,off,qoff)    

C (original Comment)
c     rotate top 'half' of cluster about z axis
c     by random angle in the range 2*pi

      theta = ran3(idum)*2*pi	! Pick random angle
      
      rot11 = cos(theta)	! Set functions of this variable
      rot12 = sin(theta)	!
      rot21 = -sin(theta)	!
      rot22 = cos(theta)	!

      do i=1,natoms/2
C For half the atoms rotate around z axis
         x = off(i)*rot11 + off(i+natoms)*rot12 ! Rotate x coordinate
         y = off(i)*rot21 + off(i+natoms)*rot22	! Rotate y coordinate
         off(i) = x 		! Set new values
         off(i+natoms) = y	! Set new values
      end do

      return	! Return
      end

C Subroutine Mutate Replace
c     Mutates one of the new strings made during mating to a new cluster

      SUBROUTINE mutate_replace(natoms,na,nb,qa,qb,r_ab,off,qoff,idum,
     &           range,nels)
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer na	! Number of atoms a 
      integer nb	! Number of atoms b
      integer idum	! Seed for random number
      integer nels	! Number of elements
C Internal integer values
      integer i
      integer iatom
C Inherited real values
      real*8 qa
      real*8 qb
      real*8 r_ab
      real *8 off(3*natoms)
      real*8 qoff(natoms)
      real*8 range
C Internal real numbers
      real*8 ran3
C Import random number generator
      external ran3

c     Generates new cluster
     
      do i=1,3*natoms
C For all atoms coordinates
         off(i) = r_ab*ran3(idum)*range
C Set coordinate as bimetallic bond length * random number * range
      end do

      if ( nels .ge. 2 ) then
C If bimetallic

         do i=1,3*natoms
C For all atoms coordinates
            off(i) = r_ab*ran3(idum)*range 
C Set coordinate as bimetallic bond length * random number * range
         end do

         do i=1,na
C Set atoms of type a
            qoff(i) = qa
         end do
         
         do i=na+1,na+nb
C Set atoms of type B
            qoff(i) = qb
         end do

      end if

      return	! Return
      end
 
C Subroutine mutate_exchange
C mutates cluster by swapping atom types
      subroutine mutate_exchange(natoms,na,nb,off,qoff,idum,pfunc,
     &                           meswaps)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer na	! Number of atoms type A
      integer nb	! Number of atoms type B
      integer idum	! Random number seed
      integer pfunc	!
      integer meswaps	! Number of mutations to perform
C Internal integer values
      integer i
      integer ncat
      integer nan
      integer cat(nmax)
      integer an(nmax)
      integer tq
C Inherited real values
      real*8 off(3*natoms)
      real*8 qoff(natoms)
C Internal Real values
      real*8 ran3
C Import random number generator
      external ran3

C (Original Comment)
c     create list of atoms types

      ncat = 0
      nan = 0

      if ( pfunc .eq. 3 .or. pfunc .eq. 5 ) then
c     (Original comment) rigid-ion potential

          do i=1,natoms
C For all atoms
             if (qoff(i) .gt. 0.d0) then
C If charge is > 0 is type A
                cat(ncat+1) = i
                ncat = ncat + 1
             else
C Else type B
                an(nan+1) = i
                nan = nan+1
             end if
          end do
   
      else if ( pfunc .eq. 4 .or. pfunc .eq. 6 ) then
c     gupta potential
         
          do i=1,natoms
c     both `charges' are positive: +1 and +2
             if (qoff(i) .eq. 1.d0) then
C If charge value = 1 it is type A
                cat(ncat+1) = i
                ncat = ncat + 1
             else
C Else is type B
                an(nan+1) = i
                nan = nan+1
             end if
          end do

      end if

      do i=1,meswaps
C For meswaps atoms

C Debug
C         write(*,'(''Swapping atoms'')')

C (Original Comment)
c     pick cation and anion

         ncat = ran3(idum)*na+1	! Cation
         nan = ran3(idum)*nb+1	! Anion

c     swap types

         tq = qoff(cat(ncat))	! Swap charges
         qoff(cat(ncat)) = qoff(an(nan))
         qoff(an(nan)) = tq
      
         cat(ncat) = nan 	! Swap types in cluster
         an(nan) = ncat

      end do

      return	! Return
      end
      
C Subroutine mutate_project
c     mutates cluster by projecting atoms on to the surface of a sphere
c     of radius n**1/3

      subroutine mutate_project(natoms,off,pi,range)
      implicit none

C Inherited integer values
      integer natoms
C Internal integer values
      integer i
C Inherited real values
      real*8 off(3*natoms)
      real*8 pi
      real*8 range
C Internal real values
      real*8 r
      real*8 sumx
      real*8 sumy
      real*8 sumz
      real*8 theta
      real*8 phi
      real*8 phi1

C (Original Comments)
c     scale coordinates to centre of mass
      sumx = 0.d0
      sumy = 0.d0
      sumz = 0.d0

      do i=1,natoms
C For all atoms, sum coordinate values
         sumx = sumx + off(i) 
         sumy = sumy + off(i+natoms)
         sumz = sumz + off(i+2*natoms)
      end do

      sumx = sumx/dble(natoms)	! Average of x coordinate
      sumy = sumy/dble(natoms)	! Average of y coordinate
      sumz = sumz/dble(natoms)	! Average of z coordinate

      do i=1,natoms
C For all atoms, move to centre
         off(i) = off(i) - sumx
         off(i+natoms) = off(i+natoms) - sumy
         off(i+2*natoms) = off(i+2*natoms) -sumz
      end do

      do i=1,natoms
C For all atoms
C (Original Comment)
c     calculate r - distance from centre of mass

         r = dsqrt(off(i)*off(i) 
     &        + off(i+natoms)*off(i+natoms) 
     &        + off(i+2*natoms)*off(i+2*natoms))

C (Original Comments)
c     calculate theta and phi
c     two possible values of phi due to phase difference

         theta = dacos(off(i+2*natoms)/r)		! ArcCos(z/r)
         phi = dacos(off(i)/(r*dsin(theta)))		! ArcCos(x/r.sin(ArcCos(z/r)))
         phi1 = dasin(off(i+natoms)/(r*dsin(theta)))	! ArcSin(y/r.sin(ArcCos(z/r)))

C (Original Comments)
c     calculate new coordinates

         off(i) = range*dsin(theta)*dcos(phi)		! x coordinate
         off(i+natoms) = range*dsin(theta)*dsin(phi1)	! y coordinate
         off(i+2*natoms) = range*dcos(theta)		! z coordinate

      end do

      return	! Return
      end
