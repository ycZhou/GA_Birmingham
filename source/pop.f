c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C (Original Comments)
c     routines that operate on population:
c     esort - sorts population in order of increasing energy
c     maintain - removes clusters of similar energy
c     term - ends ga run if population remains the same for n generations
c     highE_mutant - adds mutants to next gen regardless of their energy
c     e_pred - eats cluters lower than an input energy
c     nextNeigh - calculates the Coordination Number of each atom in each cluster


      SUBROUTINE esort(natoms,nstride,pop,q,energy,mutant,length)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none
C Inherited Integer values
      integer natoms	! Number of atoms
      integer nstride
      integer length	!
      integer mutant(length)	! Mutants
C Internal Integer values
      integer nels	! Number of elements
      integer ir	!
      integer i		!
      integer j		!
      integer k		!
      integer l		!
      integer tmut	!
C Inherited Real Values
      real*8 pop((3*natoms+1)*length)
      real*8 q(natoms*length)
      real*8 energy(length)
C Internal Real Values
      real*8 tarr1		! Temp Energy
      real*8 tclust(3*nmax+1)	! Temp cluster coordinates
      real*8 tq(nmax)		! Temp charge

C (Original Comment)
c     Sorts the energy array in order of increaing energy
c     and makes the corresponding re-arrangment of other arrays
c     Uses the Numerical Recipes heap sort

c     Start of heap sort

      l = (length)/2+1		! Middle point
      ir = length		! Duplicated
      
 10   continue
      
      if(l.gt.1)then			! if array size is greater than 1

         l=l-1				! Take one below middle point

         tarr1=energy(l)		! Get temp parent energy
         tmut=mutant(l)			! Get temp mutant energy

         do k=1,3*natoms+1
C For all atoms copy cluster coordinates (x,y,z) to temp array
            tclust(k) = pop(k+(l-1)*nstride)
         end do

         do k=1,natoms
C For all atoms copy charges across to temp array
            tq(k) = q(k+(l-1)*natoms)
         end do

      else				! else if array size is smaller than 1

         tarr1=energy(ir)		! Get temp parent energy
         tmut=mutant(ir)		! Get temp mutant energy

         do k=1,3*natoms+1
C For all atoms copy cluster coordinates (x,y,z) to temp array
            tclust(k)=pop(k+(ir-1)*nstride)
         end do

         do k=1,natoms
C For all atoms copy charges across to temp array
            tq(k) = q(k+(ir-1)*natoms)
         end do

         energy(ir)=energy(1)		! Set energy value at end to first
         mutant(ir)=mutant(1)		! Set mutant energy value as end to first

         do k=1,3*natoms+1
C For all atoms copy cluster coordinates (x,y,z) to new position
            pop(k+(ir-1)*nstride) = pop(k)
         end do

         do k=1,natoms
C For all atoms copy charges across to temp array
            q(k+(ir-1)*natoms) = q(k)
         end do
C For all atoms copy charges across to temp array
         ir=ir-1		! Subtract 1 from ir
         if(ir.eq.1)then	! Check if we are at last value

            energy(1)=tarr1	! Set Energy to temp value
            mutant(1)=tmut	! Set mutant energy to temp value

            do k=1,3*natoms+1
C For all atoms copy cluster coordinates (x,y,z) to new position
               pop(k) = tclust(k)
            end do

            do k=1,natoms
C For all atoms copy charges across to temp array
               q(k) = tq(k)
            end do

            return		! Return function
         endif
      endif

      i=l		! Set i to middle-1
      j=l+l		! Set j to middle

 20   if(j.le.ir)then			! if middle <= length
         if(j.lt.ir)then		! if middle < length
            if(energy(j).lt.energy(j+1))j=j+1		! swap j if energy of next is greater
         endif

         if(tarr1.lt.energy(j))then
C If temp energy is greater than that at middle(+1)

            energy(i)=energy(j)		! Swap energies
            mutant(i)=mutant(j)		! Swap energies
            do k=1,3*natoms+1
C Swap coordinated
               pop(k+(i-1)*nstride) = pop(k+(j-1)*nstride)
            end do

            do k=1,natoms
C Swap charges
               q(k+(i-1)*natoms) = q(k+(j-1)*natoms)
            end do

            i=j		! middle(+1)
            j=j+j	! Double j
         else
C Otherwise j = length + 1
            j=ir+1
         endif
C Return to start of loop
         go to 20
      endif

      energy(i)=tarr1		! Energy is saved from temp
      mutant(i)=tmut		! Energy is saved from temp

      do k=1,3*natoms+1
C Swap coordinate
         pop(k+(i-1)*nstride) = tclust(k)
      end do

      do k=1,natoms
C Swap charges 
         q(k+(i-1)*natoms) = tq(k)
      end do

C Return to start of subroutine
      go to 10

      END

C
C Subroutine maintain
C Removes clusters of similar energies
C
      
      SUBROUTINE maintain(natoms,nclust,nstride,nels,pop,q,energy,
     &     length,pfunc,GAconv)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer nclust	! Number of clusters
      integer nstride	! Step size across coordinates
      integer nels	! Number of elements
      integer pfunc	! Potential type
      integer length	! Length of array
C Internal integer values
      integer i
      integer j
      integer mark(ncmax)	! Temporary array to store deleted positions
      integer count		! Temporary counter
      integer diff              ! difference between nclust and count
      integer diffcount         ! Temporary counter
      integer temp              ! Temporary storage
C Inherited real values
      real*8 pop((3*natoms+1)*length)	! Population
      real*8 q(natoms*length)		! Atomic charge
      real*8 energy(length)		! Energies
      real*8 Edel              ! Energy range identifying similar structures
      real*8 GAconv            ! Convergence criteria of the GA

c     Removes duplicate strings from the population

      do i=1,length
C For all positions set reference to 0
         mark(i) = 0
      end do
      
      count = 1		! Counter set to start
      mark(1) = 1	! First point in array set to start

! Sven I've moved this outside the loop.
! Doesn't need to be set every time round as Edel doesn't change
c for DFT runs
      if (pfunc .eq. 6) then
         Edel = GAconv * 2.d0
c if not DFT calculation
      else
         Edel = 1.e-5
      end if
c added by SH 03/11 for DFT/GA runs

      do i=2,length
C For entire length of array beyond first atom
C (Original Comment)         
c     Identify unique strings

         if((energy(i)-energy(i-1)).gt. Edel) then
C If energies are essentially the same as previous point
C Note: this is used after arrays have been sorted
            count = count + 1		! Increase count of identical strings by 1
            mark(count) = i		! Set next duplicate position to counter point

C Debug
c            print *,'E= ',energy(i),' count= ',count,' i= ',i
         
         end if
            
      end do

C This will become redundant

      if ( count .lt. nclust) then
C Check if number of clusters is less than counter
C Print error if this (oddly) happens
C This would only be possible if length > nclust

         write(*,
     &      '(''Warning: too many clusters removed by maintain'',/)')
C added by SH 03/11. Ensures that the generation size is conserved         
         write(*,
     &       '(''Preserving population size 
     &           by limiting number of clusters removed'',/)')

C How many cluster are missing in the generation?

         diff = nclust - count
         diffcount = 0

C search for the lowest lying similar clusters

         do i =2,length
           if((energy(i)-energy(i-1)).le. Edel) then
             mark(count+i-1) = i
             diffcount = diffcount + 1
           end if
           if (diffcount .eq. diff) exit
         end do

C sort the mark values

         do i = 1,nclust-1
           do j = i, nclust
             if (mark(i) .gt. mark(j)) then
               temp = mark(i)
               mark(i) = mark(j)
               mark(j) = temp
             end if
           end do
         end do

C set the count value to it minimum size
         
         count = nclust

      end if 

C (Original Comment)
c     Remove non-unique strings

      do i = mark(1),count
C For point need to be removed to end of counter
C Loop through all possible options and move them forward by one

         energy(i) = energy(mark(i))	! Energy is energy from next individual point
         do j = 1,3*natoms+1
C Move forward coordinates
            pop(j+(i-1)*nstride) = pop(j+(mark(i)-1)*nstride)
         end do

         if ( nels .ge. 2 ) then
C If number of elements is greater than 2
            do j = 1,natoms
C Move charges forward
               q(j+(i-1)*natoms) = q(j+(mark(i)-1)*natoms)
            end do
         end if
      end do

      return 	! End of subroutine
      end

C
C Subroutine Term
C ends ga run if population remains the same for n generations
C

      subroutine term(nclust,iterm,tcount,energy,tflag,GAconv,
     &                oldcluster)
      implicit none

C Inherited integer values
      integer nclust	! Number of clusters
      integer iterm	! Iteration Count with this value
      integer tcount	! Termination Count
C Internal integer values
      integer i
C Inherited real values
      real*8 energy(nclust)	! Energy of clusters (Why do we need this?)
C Internal real values
      real*8 rcurr		! Current Energy Difference to preious min
      real*8 rprev		! Previous Energy Min
      real*8 GAconv             ! GA convergence criteria
      real*8 oldCluster         ! Lowest energy in old generation
C Termination Bollean flag
      logical tflag		! Terminate?
      
C Save previous value in function
      save rprev
C Set to 0
      data rprev /0.d0/

C (Original Comment)
c     test for convergence
      
      rcurr = abs(energy(nclust) - energy(1))	! Subtract energy of end cluster from first to see if they are the same
! If not DFT normal termination routine
      if (GAconv .eq. 1.d-6) then
        ! GA Emperical
        rcurr = abs(energy(nclust) - energy(1)) ! Subtract energy of end cluster from first to see if they are the same
      else
        rcurr = abs(oldCluster - energy(1)) ! Compare generation for DFT
! This if statement is unecessary as the minimum energy will never increase. We'll just set it each time
C       if (energy(1) .lt. oldCluster) oldCluster = energy(1)
        oldCluster = energy(1)
      end if
      ! This is if the whole generation is the same
      ! Otherwise rcurr is the breadth of energies in the population
      if (rcurr .lt. GAconv) rcurr = 0.d0       ! If so set residual to 0

      if (rcurr .eq. rprev) then
C If residual of cluster comparison is 0
         iterm = iterm + 1              ! Increment iteration
      else
C Else set to zero
       iterm = 0
      end if
     
C      print *, 'iterm=',iterm

C Set previous energy difference to current
      rprev = rcurr

C (Original comment)
c     stop if no improvement for tcount generations

      if (iterm .eq. tcount) tflag = .true.     ! Terminate

      return    ! End

      end 

C
C Subroutine highE_mutant
C adds mutants to next gen regardless of their energy
C

      subroutine highE_mutant(natoms,nclust,nstride,noff,idum,nels,
     &     pop,q,mutant,energy)
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms in each cluster
      integer nclust	! Number of clusters
      integer nstride	! Step size
      integer noff	! Number of Offspring
      integer idum	! Random number seed
      integer nels
      integer mutant(nclust+noff)
C Internal C values
      integer i
      integer fmut
C Inherited real values
      real*8 pop((3*natoms+1)*(nclust+noff))
      real*8 q(natoms*(nclust+noff))
      real*8 energy(nclust+noff)
C Internal real (Random number)
      real*8 ran3
C Internal boolean
      logical mflag
C Import random number generator
      external ran3

      if (ran3(idum) .lt. 0.1d0) then
C If random number is less than 0.1 (i.e. 10% chance)
         fmut = 0		! Initialise
	 mflag = .false.	! Mutant flag set false
         
         do i=nclust+1,nclust+noff
C For number of mutant clusters
            if (mutant(i) .eq. 2 .and. .not. mflag) then
C If mutant value is 2 and not mutated yet
               fmut = i		! Mark as cluster to mutate
               mflag = .true.	! Set as identifier as true
            end if
         end do

         if (mflag) then
C If we have a mutant selected
            write(*,'(''Adding higher energy mutant'',/)')

            do i=1,3*natoms+1
C For all atoms from last cluster move forward coordinates to remove cluster to mutate
               pop(i+(nclust-1)*nstride) = pop(i+(fmut-1)*nstride)
            end do
            
            if ( nels .ge. 2 ) then
C If bimetallic...
               do i=1,natoms
C Move charges also...
                  q(i+(nclust-1)*natoms) = q(i+(fmut-1)*natoms)
               end do
            end if

C Move energy values
            energy(nclust) = energy(fmut)	! Copy energy from mutant to last value
            mutant(nclust) = mutant(fmut)	! Copy mutant value to last value

C Debug
c	    write(*,'(''Energy = '',f11.6)')energy(nclust)

         end if

      end if

      return	! Return
      end

C
C Subroutine e_pred
C High Energy predator (i.e. remove bad clusters)
C

      subroutine e_pred(natoms,nclust,nstride,noff,nels,pop,q,
     &     energy,gmin)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer nclust	! Number of clusters
      integer nstride	! Step size (atoms*3)
      integer noff	! Number of offspring
      integer nels	! Number of elements
C Internal integer values
      integer i
      integer j
      integer ilive
      integer mark(ncmax)
C Inherited real values
      real*8 pop((3*natoms+1)*(nclust+noff))
      real*8 q(natoms*(nclust+noff))
      real*8 energy(nclust+noff),gmin

      do i=1,ncmax
C For all atoms, set count point to zero
         mark(i) = 0
      end do

      ilive=0	! Set reference to zero

C (Original comment)
c     find clusters with lower energy

      do i=1,nclust+noff
C For all clusters

         if (energy(i) .gt. (gmin+1.d-6)) then
C If energy is greater than the gradient minimum (i.e. not converged)

            ilive = ilive+1	! Counter ++
            mark(i) = ilive	! Add counter to array of marks
         end if

C Debug
c         print *,ilive,mark(i),energy(i),gmin

      end do
         
      if (ilive .lt. nclust) then
C If counter is less than number of clusters
         pause 'Extinction.'
      end if
      
C (Original Comment)
c     remove these clusters

      do i=1,nclust+noff
C For all clusters in offspring and parent population

         if (mark(i) .ne. 0) then
C If array of is marked with counter value

            do j=1,3*natoms+1
C For all atoms move forward one
               pop(j+(mark(i)-1)*nstride) = pop(j+(i-1)*nstride)
            end do

            if (nels .ge. 2) then
C If bimetallic
               do j=1,natoms
C Move charges forward one
                  q(j+(mark(i)-1)*natoms) = q(j+(i-1)*natoms)
               end do
            end if

            energy(mark(i)) = energy(i)		! Move energy to correct position
         end if
      end do 

      return	! Return
      end

CCCCCCCCCCCCCCCC REDUNDANT CCCCCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE nextNeigh(natoms,nclust,nstride,nels,pop,q,
C     &     length,r_ab,i)
C      implicit none
C      include "ga.param"
         
C Inherited integer values
C      integer natoms    ! Number of atoms
C      integer nclust    ! Number of clusters
C      integer nstride   ! Step size across coordinates
C      integer nels      ! Number of elements
C      integer pfunc     ! Potential type
C      integer length    ! Length of array
C Internal integer values
C      integer i         ! Generation counter from ga.f
C      integer j
C      integer k
C      integer l
C      integer mark(ncmax)       ! Temporary array to store deleted positions
C      integer count             ! Temporary counter
C      integer diff              ! difference between nclust and count
C      integer diffcount         ! Temporary counter
C      integer temp              ! Temporary storage
C Inherited real values
C      real*8 pop((3*natoms+1)*length)   ! Population
C      real*8 q(natoms*length)           ! Atomic charge
C      real*8 dist                ! Distance to next neighbour
C      real*8 nneigh(nclust,natoms) ! Storing the number of next neighbours
C      real*8 r_ab                ! Bond length as defined in input
C Character values converted
C      character*5 i_char         ! Converted i counter
C      character*50 filename      
C Logical values
C      logical fexist

C      write(i_char,'(I4)')i
C      nneigh = 0.d0
C      dist = 0.d0  

C      do k = 1,nclust
C        do j = 1,natoms
C          do l = 1,natoms
C Calculate all possible distances in cluster
C            dist = sqrt((pop(l+nstride*(k-1))-pop(j+nstride*(k-1)))**2
C     &+(pop(l+natoms+nstride*(k-1))-pop(j+natoms+nstride*(k-1)))**2
C     &+(pop(l+2*natoms+nstride*(k-1))-pop(j+2*natoms+nstride*(k-1)))**2)
C Is atom a nearest neighbour?
C            if (dist .le. r_ab) nneigh(k,l) = nneigh(k,l) + 1.d0 
C          end do
C        end do
C Substract the 1 because it counted the self term
C        nneigh = nneigh - 1.d0
C Sort numbers. Could slow down code for large clusters. No
C problem for small clusters.
C        do l = 1,natoms
C          do j = 1,natoms
C            if (nneigh(k,l) .gt. nneigh(k,j)) then
C              temp = nneigh(k,l)
C              nneigh(k,l) = nneigh(k,j)
C              nneigh(k,j) = temp
C            end if
C          end do
C        end do
C If data is only written to file this does not have to be
C used yet.
C      end do  
      
C      filename='Generation_'//adjustl(trim(i_char))
       
C      inquire(file=trim(filename),exist=fexist)
C       
C      if (fexist) then
C        open(unit=99,file=trim(filename),status='old')
C        close(unit=99,status='delete')
C      else
C        open(unit=99,file=trim(filename),status='new')
C        do j = 1,nclust
C        write(99,'(I4,2X,F12.8)')j,pop(nstride*j)
C          do k = 1,natoms
C            write(99,'(F5.3,2X,F12.8,2X,F12.8,2X,F12.8,2X,F5.3)')
C     &      q(k+natoms*(j-1)),
C     &      pop(k+nstride*(j-1)),pop(k+natoms+nstride*(j-1)),
C     &      pop(k+2*natoms+nstride*(j-1)),nneigh(j,k)
C          end do
C        end do
C        close(unit=99,status='keep')
C      end if

C      END SUBROUTINE nextNeigh
