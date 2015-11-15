c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

c     routines for selectin parents
c     roul-select - roulette wheel
c     tourn-select - tournament

C Subroutine roul_select

      SUBROUTINE roul_select(natoms,nclust,nstride,nels,pop,q,fit,pair,
     &             qp,idum,str)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer nclust	! Number of clusters
      integer nstride	! Number of steps between each cluster
      integer nels	! Number of elements
      integer idum	! Random seed
      integer str(2)
C Internal integer values
      integer i
      integer istring
      integer isel
      integer mark(ncmax)
      integer iarr
C Inherited real values
      real*8 pop((3*natoms+1)*nclust)
      real*8 q(natoms*nclust)
      real*8 fit(nclust)
      real*8 pair(6*natoms)
      real*8 qp(2*natoms)
C Internal real values
      real*8 ran3
      real*8 rd
C Import random number
      external ran3

c     Selects are pair of strings for mating
c     from the first nclust clusters of the
c     population

c     Selects a pair of strings  based on their probability
c     and copies them into pair
            
      do i=1,nclust
C For all clusters set mark as none
         mark(i) = 0
      end do

      isel = 0	! Set selection

      do while(isel .lt. 2)
C Whilst selected is less than two

         istring = ran3(idum)*nclust+1	! Pick position for random clusters
         rd = ran3(idum)		! Random number

         if (rd.lt.fit(istring) .and. mark(istring).ne.1) then
C If fitness is greater than random number, and not already selected

            iarr = isel*natoms	! Create array for coordinates
            mark(istring) = 1	! Mark value as 1
            do i=1,3*natoms
C For all atom coordinates
C Not sure what is going on here. Selecting coordinates of current atom
               pair(i+3*iarr) = pop(i+(istring-1)*nstride)  
            end do

            if (nels .gt. 1) then
C If number of elements > 1
               do i=1,natoms
C For all atoms copy charge of selected atom across
                  qp(i+iarr) = q(i+(istring-1)*natoms)
               end do
            end if

            isel = isel + 1		! Number selected ++
            str(isel) = istring		! Turn into string
         end if
      end do

C Debug
c      write(7,
c     &'('' Str1 = '',i2,'' f = '',f6.4,'' Str2 = '',i2,'' f = '',f6.4)')
c     &str(1),fit(str(1)),str(2),fit(str(2))

      return	! Return
      end

C Subroutine tourn_select
c     Selects are pair of strings for mating
c     from the first nclust clusters of the
c     population


      SUBROUTINE tourn_select(natoms,nclust,nstride,nels,pop,q,energy,
     &            pair,qp,idum,tourn_size,str)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer nclust	! Number of clusters
      integer nstride
      integer nels
      integer idum
      integer tourn_size	! Tournament size
      integer str(2)
C Internal integer values
      integer i
      integer j
      integer best
      integer second
      integer tourn(ncmax)
      integer mark(ncmax)
      integer n
      integer iran
      integer temp
C Inherited real values
      real*8 pop((3*natoms+1)*nclust)	! Coordinates
      real*8 q(natoms*nclust)		! Charge
      real*8 energy(nclust)		! Energy
      real*8 pair(6*natoms)		! Pair of atoms
      real*8 qp(2*natoms)		! Pair charges
C Internal real values
      real*8 ran3
C Export random number generator
      external ran3

      do i=1,nclust
C For all clusters set marks to zero
         mark(i) = 0
      end do

C (Original Comment)
c     Selects a pair of strings using a tournament selection
c     pick tourn_size strings at random
      
      n = 0	! Counter of tournament size

      do while(n .lt. tourn_size)
C Whilst population not filled

         iran = ran3(idum)*nclust+1	! Get random number for cluster position
         if (mark(iran) .ne. 1) then
C If not already selected
            n = n+1		! Add one to population
            tourn(n) = iran	! Copy position into tournament
            mark(iran) = 1	! Set marked position
         end if
      end do

C (Original Comments)
c     find two lowest energy clusters
c     straight insertion sort from Numerical Recipes

      do i=2,tourn_size
C For entire tournament
         temp = tourn(i)
         do j=i-1,1,-1
C Compare current energy to others around
            if (energy(tourn(j)) .le. energy(temp)) goto 10
            tourn(j+1) = tourn(j)	! Move value backwards
         end do
         j=0				! Reset j
 10      tourn(j+1) = temp		! Move energy value in array
      end do

      best = tourn(1)		! Best value
      second = tourn(2)		! Second best value

      str(1) = best		! Change to strings ...
      str(2) = second		!

      do i=1,3*natoms
C For all atoms in cluster set value in pair array
         pair(i) = pop(i+(best-1)*nstride)
      end do
      
      do i=1,3*natoms
C Set second set of atoms as that to be paired with the first
         pair(i+nstride) = pop(i+(second-1)*nstride)
      end do

      if (nels .gt. 1) then
C If we have a bimetallic system

         do i=1,natoms
C For all atoms copy charges of best cluster
            qp(i) = q(i+(best-1)*natoms)
         end do

         do i=1,natoms
C For all atoms copy charges of second best cluster
            qp(i+natoms) = q(i+(second-1)*natoms)
         end do

      end if

C Debug
c      write(*,'(a,1x,a,i2,a,1x,a,i2,a)')'sel','a',str(1),'a',
c     &     'a',str(2),'a'
      
      return	! Return
      end
      
      
