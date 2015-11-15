c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     (Original comments)
c     routines for generating initial population :
c     gen_pop_box - box of N^(1/3)
c     gen_pop_sphere  - puts atoms of surface of sphere of radius N^(1/3)

      subroutine gen_pop_box(natoms,nclust,nstride,pop,idum,range)
      implicit none

C Integer Values
      integer natoms			! Number of atoms
      integer nclust			! Number of clusters
      integer nstride			! 
      integer idum			! Random Number Seed
      integer i
      integer j
C Real Values
      real*8 pop((3*natoms+1)*nclust)	! Population
      real*8 range
      real*8 ran3
C Import random number generator
      external ran3

C (Original comments)
c     seed random number generator and generate initial population

C Loop through clusters
      do i=1,nclust
C Loop through atom coordinates
         do j=1,3*natoms
            pop(j+(i-1)*nstride)=ran3(idum)*range
         end do

      end do

      return
      end

C Begin subroutine for cluster generation on sphere
      subroutine gen_pop_sphere(natoms,nclust,nstride,pop,idum,range,pi)
      implicit none

C Integer Values
      integer natoms                    ! Number of atoms
      integer nclust                    ! Number of clusters
      integer nstride                   ! 
      integer idum                      ! Random Number Seed
      integer i
      integer j
C Real Values
      real*8 pop((3*natoms+1)*nclust)
      real*8 range
      real*8 pi
      real*8 ran3
      real*8 srange
      real*8 theta
      real*8 phi
      real*8 dist2
C Import random number generator
      external ran3

C (Original Comments)
c     Generate random positions on the surface of a sphere
c     set radius of sphere to (n**1/3)/2 - produces cages on minimisation

      srange = range/2.d0
      
C Loop over clusters
      do i=1,nclust
C Loop over atom x coordinates 
         do j=1,natoms
            theta = ran3(idum)*2.d0*pi
            phi = ran3(idum)*2.d0*pi
            pop(j+(i-1)*nstride) = srange*dsin(theta)*dcos(phi)
         end do
C Loop over atom y coordinates
         do j=natoms+1,2*natoms
            theta = ran3(idum)*2.d0*pi
            phi = ran3(idum)*2.d0*pi
            pop(j+(i-1)*nstride) = srange*dsin(theta)*dsin(phi)
         end do
C Loop over atoms z coordinates            
         do j=2*natoms+1,3*natoms
            theta = ran3(idum)*2.d0*pi
            phi = ran3(idum)*2.d0*pi
            pop(j+(i-1)*nstride) = srange*dcos(theta)
         end do

      end do

      return
      end
