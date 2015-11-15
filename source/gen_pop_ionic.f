c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C Subroutine to generate population for bimetallic clusters
      subroutine gen_pop_ionic(na,nb,nclust,nstride,idum,pop,q,qa,qb,
     &           r_ab,range)
      implicit none

C Integer Values
      integer na
      integer nb
      integer nclust		! Number of Clusters
      integer nstride
      integer idum		! Random Number Generator Seed
      integer i
      integer j
C Real Values
      real*8 pop((3*(na+nb)+1)*nclust)	! Population
      real*8 q((na+nb)*nclust)		! Array of charges
      real*8 qa				! Charge of a
      real*8 qb				! Charge of b
      real*8 r_ab			! Distance between a and b
      real*8 range			! Range of atoms 
      real*8 ran3			! Random number

      external ran3

C     (Original comment)
c     seed random number generator and generate initial population
C Cycle through number of clusters
      do i=1,nclust
C Loop through all atoms and assign coordinates         
         do j=1,3*(na+nb)
            pop(j+(i-1)*nstride) = r_ab*ran3(idum)*range
         end do
Calculate charge for atom type A
         do j=1,na
            q(j+(i-1)*(na+nb)) = qa
         end do
Calculate charge for atom type B
         do j=na+1,na+nb
            q(j+(i-1)*(na+nb)) = qb
         end do
      end do

      return
      end
