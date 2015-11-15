c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

c     mating routines :
c     mate_1pt_random - 1pt crossover about random postion
c     mate_1pt_weighted - 1pt crossover about position derived from relative fitness values
c     mate_2pt
c     mate_inter - averages coordinates of parents
c     mate_2pt_multi_random
c     mate_exchange - swaps atom types in ionic clusters
 
C Subroutine mate_1pt_random
c     Cuts and pastes two selected strings about a single random point

      SUBROUTINE mate_1pt_random(natoms,off,pair,idum)
      implicit none

C Inherited integer variables
      integer natoms	! Number of atoms
      integer idum	! Random number seed
C Internal integer variables
      integer i
      integer posit
      integer nvar
C Inherited real number variables
      real*8 off(3*natoms)
      real*8 pair(6*natoms)
C Internal real variable            
      real*8 ran3
C Import random number generator
      external ran3

      nvar = 3*natoms			! Number of coordinates
      posit = ran3(idum)*(natoms-1) + 1	! Random times position 

      do i=1,posit
         off(i) = pair(i)
         off(natoms+i) = pair(natoms+i)
         off(2*natoms+i) = pair(2*natoms+i)
      end do

      do i=posit+1,natoms
         off(i) = pair(i+nvar)
         off(natoms+i) = pair(natoms+i+nvar)
         off(2*natoms+i) = pair(2*natoms+i+nvar)
      end do
       
C So we end up with linear combination of A and B

      return	! Return
      end

C Subroutine mate_1pt_weighted

      subroutine mate_1pt_weighted(natoms,off,pair,fit1,fit2)
      implicit none

C Inherited integer values
      integer natoms
C Internal integer values
      integer i
      integer posit
      integer nvar
C Inherited real values
      real*8 off(3*natoms)
      real*8 pair(6*natoms)
      real*8 fit1	! Fitness Cluster 1
      real*8 fit2	! Fitness Cluster 2
      
      nvar = 3*natoms			! Number of coordinates
      posit = natoms*(fit1/(fit1+fit2))	! Weighting with respect to fitness of clusters

C Debug
c      print *,posit,fit(str(1)),fit(str(2))

C As previously, copy parts of clusters as required
      do i=1,posit
         off(i) = pair(i)
         off(natoms+i) = pair(natoms+i)
         off(2*natoms+i) = pair(2*natoms+i)
      end do

      do i=posit+1,natoms
         off(i) = pair(i+nvar)
         off(natoms+i) = pair(natoms+i+nvar)
         off(2*natoms+i) = pair(2*natoms+i+nvar)
      end do
       
      return	! Return
      end

C Subroutine mate_2pt
Cc     performs 2 point crossover about random positions
      
      SUBROUTINE mate_2pt(natoms,off,pair,idum)
      implicit none

      integer natoms,idum
      real*8 off(3*natoms),pair(6*natoms)

      integer i,pos1,pos2,nvar
      real*8 ran3

      external ran3

C (Original Comments)
c     pick positions
      
      nvar = 3*natoms				! Total Coordinates
      pos1 = ran3(idum)*(natoms-2)+1		! Random Position 1
      pos2 = ran3(idum)*(natoms-(pos1+1)) + pos1 + 1	! Random Position 2
      
C (Original Comment)
c     cut and paste about positions

      do i=1,pos1
         off(i) = pair(i)
         off(natoms+i) = pair(natoms+i)
         off(2*natoms+i) = pair(2*natoms+i)
      end do

      do i=pos1+1,pos2
         off(i) = pair(i+nvar)
         off(natoms+i) = pair(natoms+i+nvar)
         off(2*natoms+i) = pair(2*natoms+i+nvar)
      end do

      do i=pos2+1,natoms
         off(i) = pair(i)
         off(natoms+i) = pair(natoms+i)
         off(2*natoms+i) = pair(2*natoms+i)
      end do

      return
      end

C Subroutine mate_inter
c     Mates by interpolating the two clusters

      SUBROUTINE mate_inter(natoms,off,pair)
      implicit none

C Inherited integer values
      integer natoms
C Internal integer variables
      integer i
C Inherited real variables
      real*8 off(3*natoms)
      real*8 pair(6*natoms)

      do i=1,natoms
C For all atoms, take average of their values and save
         off(i) = (pair(i)+pair(i+3*natoms))*0.5d0
         off(i+natoms) = (pair(i+natoms)+pair(i+4*natoms))*0.5d0
         off(i+2*natoms) = (pair(i+2*natoms)+pair(i+5*natoms))*0.5d0
      end do

      return ! Return
      end

C Subroutine mate_1pt_multi_random
c     single point crossover about a random position

      SUBROUTINE mate_1pt_multi_random(natoms,pfunc,off,qoff,pair,qp,
     &     idum)
      implicit none
      
C Inherited integer variables
      integer natoms	! Number of atoms
      integer pfunc	! Function type
      integer idum	! Random Number Seed
C Internal integer values
      integer i
      integer j
C      integer ncount
      integer posit
      integer ca
      integer caa
      integer cb
      integer cbb
C      integer ta
C      integer tb
C Inherited real variables
      real*8 off(3*natoms)	! Coordinates
      real*8 qoff(natoms)	! Charges
      real*8 pair(6*natoms)	! Clusters to mate
      real*8 qp(2*natoms)
C Internal Real Variable
      real*8 ran3
C Import random number generator
      external ran3

c     generated random atom position
      posit = ran3(idum)*natoms+1

      ca = 0	! Counter type A
      cb = 0	! Counter type B

      if (pfunc .eq. 3 .or. pfunc .eq. 5) then
C If we are dealing with an ionic potential...

         do i=1,posit
C (Original Comments)
c     copy i to posit atoms from first cluster
            qoff(i) = qp(i)
c     count number of each atom type copied
            if (qoff(i) .gt. 0.d0) then
               ca = ca+1	! Type A
            else
               cb = cb+1	! Type B
            end if

C Copy atoms across
            off(i) = pair(i)
            off(i+natoms) = pair(i+natoms)
            off(i+2*natoms) = pair(i+2*natoms)
         end do
         
         caa = 0		! Second counter of atoms a
         cbb = 0		! Second counter of atoms b
         posit = posit+1	! Position is moved forward 1
      
C (Original Comments)
c     copy remaining atoms of each type from second cluster
c     ignores first ca and cb atoms in each cluster
         
         do i=1,natoms
C For all atoms
            if(qp(i+natoms) .gt. 0.d0) then
C If type A
               if (caa .lt. ca) then
C And counter is not yet equal to those found for cluster A keeping filling spaces
                  caa = caa + 1
               else
C Otherwise copy details all across
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit + 1	! Position ++
               end if
            else
C Otherwise type B
               if (cbb .lt. cb) then
C If not reached number of type B from cluster A, keep filling spaces
                  cbb = cbb + 1		! ++
               else
C Other wise transfer atom across with charge
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit+1
               end if
            end if
         end do

      else
C If not dealing with an ionic potential
         
         do i=1,posit
C (Original comment)
c     copy i to posit atoms from first cluster
            qoff(i) = qp(i)
c     count number of each atom type copied
            if (qoff(i) .eq. 1.0d0) then
               ca = ca+1
            else
               cb = cb+1
            end if
C Copy atom details across
            off(i) = pair(i)
            off(i+natoms) = pair(i+natoms)
            off(i+2*natoms) = pair(i+2*natoms)
         end do
         
         caa = 0		! Initiate second counters
         cbb = 0		!
         posit = posit+1	!
      
         do i=1,natoms
c     copy remaining atoms of each type from second cluster
c     ignores first ca and cb atoms in each cluster

            if(qp(i+natoms) .eq. 1.0d0) then
C If type A
               if (caa .lt. ca) then
                  caa = caa + 1
               else
C Copy positions and charges
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit + 1
               end if
            else
C Otherwise type B
               if (cbb .lt. cb) then
                  cbb = cbb + 1
               else
C Copy positions and charges
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit+1
               end if
            end if
         end do
         
      end if
      
      return	! Return
      end

C Subroutine mate_1pt_multi_weighted
C single point crossover about a weighted position

      SUBROUTINE mate_1pt_multi_weighted(natoms,pfunc,off,qoff,pair,qp,
     &            fit1,fit2)
      implicit none

C Inherited Integer Variables
      integer natoms
      integer pfunc
C Internal integer values
      integer i
      integer j
C      integer ncount
      integer posit
      integer ca
      integer caa
      integer cb
      integer cbb
C      integer ta
C      integer tb
C Inherited real variables
      real*8 off(3*natoms)
      real*8 qoff(natoms)
      real*8 pair(6*natoms)
      real*8 qp(2*natoms)
      real*8 fit1
      real*8 fit2

C (Original Comment)
c     generated random atomic position
      posit = natoms*(fit1/(fit1+fit2))

      ca = 0    ! Counter type A
      cb = 0    ! Counter type B

      if (pfunc .eq. 3 .or. pfunc .eq. 5) then
C If we are dealing with an ionic potential...

         do i=1,posit
C (Original Comments)
c     copy i to posit atoms from first cluster
            qoff(i) = qp(i)
c     count number of each atom type copied
            if (qoff(i) .gt. 0.d0) then
               ca = ca+1        ! Type A
            else
               cb = cb+1        ! Type B
            end if

C Copy atoms across
            off(i) = pair(i)
            off(i+natoms) = pair(i+natoms)
            off(i+2*natoms) = pair(i+2*natoms)
         end do
         
         caa = 0                ! Second counter of atoms a
         cbb = 0                ! Second counter of atoms b
         posit = posit+1        ! Position is moved forward 1
      
C (Original Comments)
c     copy remaining atoms of each type from second cluster
c     ignores first ca and cb atoms in each cluster
         
         do i=1,natoms
C For all atoms
            if(qp(i+natoms) .gt. 0.d0) then
C If type A
               if (caa .lt. ca) then
C And counter is not yet equal to those found for cluster A keeping filling spaces
                  caa = caa + 1
               else
C Otherwise copy details all across
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit + 1     ! Position ++
               end if
            else
C Otherwise type B
               if (cbb .lt. cb) then
C If not reached number of type B from cluster A, keep filling spaces
                  cbb = cbb + 1         ! ++
               else
C Other wise transfer atom across with charge
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit+1
               end if
            end if
         end do

      else
C If not dealing with an ionic potential
         
         do i=1,posit
C (Original comment)
c     copy i to posit atoms from first cluster
            qoff(i) = qp(i)
c     count number of each atom type copied
            if (qoff(i) .eq. 1.0d0) then
               ca = ca+1
            else
               cb = cb+1
            end if
C Copy atom details across
            off(i) = pair(i)
            off(i+natoms) = pair(i+natoms)
            off(i+2*natoms) = pair(i+2*natoms)
         end do
         
         caa = 0                ! Initiate second counters
         cbb = 0                !
         posit = posit+1        !
      
         do i=1,natoms
c     copy remaining atoms of each type from second cluster
c     ignores first ca and cb atoms in each cluster

            if(qp(i+natoms) .eq. 1.0d0) then
C If type A
               if (caa .lt. ca) then
                  caa = caa + 1
               else
C Copy positions and charges
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit + 1
               end if
            else
C Otherwise type B
               if (cbb .lt. cb) then
                  cbb = cbb + 1
               else
C Copy positions and charges
                  qoff(posit) = qp(i+natoms)
                  off(posit) = pair(i+3*natoms)
                  off(posit+natoms) = pair(i+4*natoms)
                  off(posit+2*natoms) = pair(i+5*natoms)
                  posit = posit+1
               end if
            end if
         end do
         
      end if
      
      return    ! Return
      end

C Update 11/11/2010 - Working fine for bimetallics
C Subroutine 2pt_multi_random
C Double point crossover about a random position

      SUBROUTINE mate_2pt_multi_random(natoms,pfunc,na,nb,off,qoff,pair,
     &            qp,idum)

      implicit none

C Inherited integer values
      integer natoms
      integer pfunc     ! Function type
      integer na
      integer nb
      integer idum
C Internal integer values
      integer i
      integer j
      integer ncount
      integer posit1
      integer posit2
      integer ca
C      integer caa	! Double check
      integer cb
C      integer cbb	! Double check
C      integer ta	! Double check
C      integer tb	! Double check
      integer empty
C Inherited real values
      real*8 off(3*natoms)
      real*8 qoff(natoms)
      real*8 pair(6*natoms)
      real*8 qp(2*natoms)
C Internal real values
      real*8 ran3
C Import random number generator
      external ran3

c     generated random atomic positions
      posit1 = ran3(idum)*(natoms-2)+1
      posit2 = ran3(idum)*(natoms-(posit1+1)) + posit1 + 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC Used to get this working CCCCCCCCCCCCCCCCCCC
C Debug    
C        do i=1,2*natoms
C Sum up total of all atoms of each type in new cluster
C           write(*,*)qp(i)
C        end do
 
C Debug
C Reset all charges and coordinates
C This is causing me problems in the 2pt crossover. Reason unknown

C          do i=1,natoms
C            qoff(i) = 0.0d0
C            off(i) = 0.0d0
C            off(i+natoms) = 0.0d0
C            off(i+2*natoms) = 0.0d0
C          end do
CC No longer needed (though some may say there is no harm in zeroing variables CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Debug
C      print *,'pos 1 ',posit1,'  pos2 ',posit2

      ca = 0	! Count of element A
      cb = 0	! Count of element B
      empty = 0	! Check for empty spaces

      if (pfunc .eq. 3 .or. pfunc .eq. 5) then
C If we are dealing with an ionic potential...

C (Original Comments)
c     copy atoms from first cluster

        do i=1,posit1
C Up to position 1
           if (qp(i) .gt. 0.d0) then
C If atom type A
              if (ca .lt. na) then
C If count of type A is less than total needed copy across
                 qoff(i-empty) = qp(i)
                 off(i-empty) = pair(i)
                 off(i-empty+natoms) = pair(i+natoms)
                 off(i-empty+2*natoms) = pair(i+2*natoms)
                 ca = ca + 1
              else
                 empty = empty + 1
              end if
           else
C Else Atom type B
              if (cb .lt. nb) then
C If count of type B is less than total needed copy across
                 qoff(i-empty) = qp(i)
                 off(i-empty) = pair(i)
                 off(i-empty+natoms) = pair(i+natoms)
                 off(i-empty+2*natoms) = pair(i+2*natoms)
                 cb = cb + 1
              else
                 empty = empty + 1
              end if
           end if
        end do

CC Redundant Code CC
C        do i=1,posit1
C Up to position 1
C           if (qp(i) .gt. 0.d0) then
C If atom type A
C              if (ca .lt. na) then
C If count of type A is less than total needed copy across
C                 qoff(i) = qp(i)
C                 off(i) = pair(i)
C                 off(i+natoms) = pair(i+natoms)
C                 off(i+2*natoms) = pair(i+2*natoms)
C                 ca = ca + 1
C              end if
C           else
C Else Atom type B
C              if (cb .lt. nb) then
C If count of type B is less than total needed copy across
C                 qoff(i) = qp(i)
C                 off(i) = pair(i)
C                 off(i+natoms) = pair(i+natoms)
C                 off(i+2*natoms) = pair(i+2*natoms)
C                 cb = cb + 1
C              end if
C           end if
C        end do
CCCCCCCCCCCCCCCCCCCCCCC

C Debug
C      print *,'ca ',ca,' cb ',cb
C (Original Comments)
c     copy atoms from second cluster

C Update posit to reflect empty spaces
        posit1 = posit1 + 1 - empty
        empty = 0

        do i=posit1,posit2
C For atoms from point 1 to point 2
           if (qp(i+natoms) .gt. 0.d0) then
C If atom Type A
              if (ca .lt. na) then
C If count of type A is less than total needed copy across
                 qoff(i-empty) = qp(i+natoms)
                 off(i-empty) = pair(i+3*natoms)
                 off(i-empty+natoms) = pair(i+4*natoms)
                 off(i-empty+2*natoms) = pair(i+5*natoms)
                 ca = ca + 1
              else
                 empty = empty + 1
              end if
           else
C Else Atom Type B
              if (cb .lt. nb) then
C If count of type B is less than total needed copy across
                 qoff(i-empty) = qp(i+natoms)
                 off(i-empty) = pair(i+3*natoms)
                 off(i-empty+natoms) = pair(i+4*natoms)
                 off(i-empty+2*natoms) = pair(i+5*natoms)
                 cb = cb + 1
              else
                 empty = empty + 1
              end if

           end if
        end do

C Debug
C      print *,'ca ',ca,' cb ',cb, ' empty ', empty

C        caa = 0         ! Count final of type A
C        cbb = 0         ! Count final of type B
        posit2 = posit2 + 1 - empty   ! Position counter

CC Redundant code CC
C        do i=posit1+1,posit2
C For atoms from point 1 to point 2
C           if (qp(i+natoms) .gt. 0.d0) then
C If atom Type A
C              if (ca .lt. na) then
C If count of type A is less than total needed copy across
C                 qoff(i) = qp(i+natoms)
C                 off(i) = pair(i+3*natoms)
C                 off(i+natoms) = pair(i+4*natoms)
C                 off(i+2*natoms) = pair(i+5*natoms)
C                 ca = ca + 1
C              end if
C           else
C Else Atom Type B
C              if (cb .lt. nb) then
C If count of type B is less than total needed copy across
C                 qoff(i) = qp(i+natoms)
C                 off(i) = pair(i+3*natoms)
C                 off(i+natoms) = pair(i+4*natoms)
C                 off(i+2*natoms) = pair(i+5*natoms)
C                 cb = cb + 1
C              end if
C           end if
C        end do

C Debug
C      print *,'ca ',ca,' cb ',cb

C        caa = 0		! Count final of type A
C        cbb = 0		! Count final of type B
C        posit2 = posit2 + 1	! Position counter
CCCCCCCCCCCCCCCCC

C Fill in remaining spaces
        i = natoms + 1

        do while (ca .lt. na .or. cb .lt. nb)
           i = i - 1
C For all atoms
           if (qp(i) .gt. 0.d0) then
C If atom type A
              if (ca .lt. na) then
C Transfer in atom details if not full
                 qoff(posit2) = qp(i)
                 off(posit2) = pair(i)
                 off(natoms+posit2) = pair(natoms+i)
                 off(2*natoms+posit2) = pair(2*natoms+i)
                 posit2 = posit2 + 1
                 ca = ca + 1
              end if
           else
C Atom type B
              if (cb .lt. nb) then
C Transfer in atom details, if we haven't peaked
                 qoff(posit2) = qp(i)
                 off(posit2) = pair(i)
                 off(natoms+posit2) = pair(natoms+i)
                 off(2*natoms+posit2) = pair(2*natoms+i)
                 posit2 = posit2 + 1
                 cb = cb + 1
              end if
           end if
        end do

CC Redundant code, succeeded by working code
C (Original Comment)
c     copy atoms from first cluster
c     atoms taken from ca + 1 and cb + 1 positions

C        do i=1,natoms
C For all atoms
C           if (qp(i) .gt. 0.d0) then
C If atom type A
C              if (caa .ge. ca) then
C If we have not reached equilibria keep adding and searching
C Else transfer in atom details
C                 if (caa .lt. na) then
C                    qoff(posit2) = qp(i)
C                    off(posit2) = pair(i)
C                    off(natoms+posit2) = pair(i+natoms)
C                    off(2*natoms+posit2) = pair(i+2*natoms)
C                    posit2 = posit2 + 1
C                 end if
C              end if
C              caa = caa + 1
C           else
C Atom type B
C              if (cbb .ge. cb) then
C If we have not reached equilibria keep adding and searching for type B
C Otherwise transfer in atom details, if we haven't peaked
C                 if (cbb .lt. nb) then
C                    qoff(posit2) = qp(i)
C                    off(posit2) = pair(i)
C                    off(natoms+posit2) = pair(i+natoms)
C                    off(2*natoms+posit2) = pair(i+2*natoms)
C                    posit2 = posit2 + 1
C                 end if
C              end if
C              cbb = cbb + 1
C           end if
C        end do

C        ta = 0		! Total Counters, atom A
C        tb = 0		! Total counters, atom B

C        do i=1,natoms
C Sum up total of all atoms of each type in new cluster
C           if (qoff(i) .gt. 0.d0) then
C              ta = ta + 1
C           else
C              tb = tb + 1
C           end if
C        end do

C        print *,'ta ',ta,' tb ',tb
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      else
C If we are dealing with a non-ionic potential... most likely gupta

C (Original Comments)
c     copy atoms from first cluster

CCCCCCCCCCCCCCC
C Debug
C        do i=1,natoms
C Sum up total of all atoms of each type in new cluster
C           write(*,*)qoff(i)
C        end do
CCCCCCCCCCCCCCC

        do i=1,posit1
C Up to position 1
           if (qp(i) .eq. 1.d0) then
C If atom type A
              if (ca .lt. na) then
C If count of type A is less than total needed copy across
                 qoff(i-empty) = qp(i)
                 off(i-empty) = pair(i)
                 off(i-empty+natoms) = pair(i+natoms)
                 off(i-empty+2*natoms) = pair(i+2*natoms)
                 ca = ca + 1
              else
                 empty = empty + 1 
              end if
           else
C Else Atom type B
              if (cb .lt. nb) then
C If count of type B is less than total needed copy across
                 qoff(i-empty) = qp(i)
                 off(i-empty) = pair(i)
                 off(i-empty+natoms) = pair(i+natoms)
                 off(i-empty+2*natoms) = pair(i+2*natoms)
                 cb = cb + 1
              else
                 empty = empty + 1
              end if
           end if
        end do

CCCCCCCCCCCCCC
C Debug
C      print *,'ca ',ca,' cb ',cb,' empty ',empty
C (Original Comments)
c     copy atoms from second cluster

C Debug
C        do i=1,natoms
C Sum up total of all atoms of each type in new cluster
C           write(*,*)qoff(i)
C        end do
CCCCCCCCCCCCCC

C Update posit to reflect empty spaces
        posit1 = posit1 + 1 - empty
        empty = 0

        do i=posit1,posit2
C For atoms from point 1 to point 2
           if (qp(i+natoms) .eq. 1.d0) then
C If atom Type A
              if (ca .lt. na) then
C If count of type A is less than total needed copy across
                 qoff(i-empty) = qp(i+natoms)
                 off(i-empty) = pair(i+3*natoms)
                 off(i-empty+natoms) = pair(i+4*natoms)
                 off(i-empty+2*natoms) = pair(i+5*natoms)
                 ca = ca + 1
              else
                 empty = empty + 1
              end if
           else
C Else Atom Type B
              if (cb .lt. nb) then
C If count of type B is less than total needed copy across
                 qoff(i-empty) = qp(i+natoms)
                 off(i-empty) = pair(i+3*natoms)
                 off(i-empty+natoms) = pair(i+4*natoms)
                 off(i-empty+2*natoms) = pair(i+5*natoms)
                 cb = cb + 1
              else
                 empty = empty + 1
              end if
              
           end if
        end do

C Debug
C      print *,'ca ',ca,' cb ',cb, ' empty ', empty

C        caa = 0         ! Count final of type A
C        cbb = 0         ! Count final of type B
        posit2 = posit2 + 1 - empty   ! Position counter

CC This is all debug from getting this part of program working CC

C        posit2 = natoms + 1
C        empty = 0

C        do i=1,natoms
C           if (qoff(i) .ne. 1.d0 .and. qoff(i) .ne. 2.d0) then
C              posit2 = posit2 - 1
C           endif
C        end do

C Debug
C        do i=1,natoms
C Sum up total of all atoms of each type in new cluster
C           write(*,*)qoff(i)
C        end do

C Debug
C      print *,'pos 1 ',posit1,'  pos2 ',posit2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C (Original Comment)
c     copy atoms from first cluster
c     start at back and fill in remaining hols

        i = natoms + 1

        do while (ca .lt. na .or. cb .lt. nb)
           i = i - 1
C For all atoms
           if (qp(i) .eq. 1.0d0) then
C If atom type A
              if (ca .lt. na) then
C Transfer in atom details if not full
                 qoff(posit2) = qp(i)
                 off(posit2) = pair(i)
                 off(natoms+posit2) = pair(natoms+i)
                 off(2*natoms+posit2) = pair(2*natoms+i)
                 posit2 = posit2 + 1
                 ca = ca + 1
              end if
C           else if (qp(i) .eq. 2.0d0) then
           else
C Atom type B
              if (cb .lt. nb) then
C Transfer in atom details, if we haven't peaked
                 qoff(posit2) = qp(i)
                 off(posit2) = pair(i)
                 off(natoms+posit2) = pair(natoms+i)
                 off(2*natoms+posit2) = pair(2*natoms+i)
                 posit2 = posit2 + 1
                 cb = cb + 1
              end if
           end if
        end do

CC This is all debug from getting this part of program working CC

C        ta = 0          ! Total Counters, atom A
C        tb = 0          ! Total counters, atom B

C        do i=1,natoms
C Sum up total of all atoms of each type in new cluster
C           if (qoff(i) .eq. 1.d0) then
C              ta = ta + 1
C           else if (qoff(i) .eq. 2.d0) then
C              tb = tb + 1
C           end if
C           write(*,*)qoff(i)
C        print *,'x ',off(i),' y ',off(i+natoms),' z ',off(i+2*natoms)
C        end do

CC
CC Reassign character types to even things out
CC This is no longer needed now problems have been resolved
CC

C Debug
C Print totals
C        print *,'ta ',ta,' tb ',tb,' na ',na,' nb ',nb
C
C       if (ta .ne. na .or. tb .ne. nb) then
C           do i=1,natoms
C           write(*,*)qoff(i)
C           print *,'x ',off(i),' y ',off(i+natoms),' z ',off(i+2*natoms)
C Sum up total of all atoms of each type in new cluster
C              if (qoff(i) .ne. 1.d0 .and. qoff(i) .ne. 2.d0) then
C                 write(*,'(''We have a problem here... No value for the
C     & element type?'')')
C                 if (ta .ne. na) then
C                    write(*,'(''Reassigning to type A'')')
C                    qoff(i) = 1.0d0
C                 else if (tb .ne. nb) then
C                    write(*,'(''Reassigning to type B'')')
C                    qoff(i) = 2.0d0               
C                 end if
C              end if
C           end do
C        end if
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end if

      return		! Return
      end

C Subroutine mate_exchange
C mates ionic clusters by exchange cations and anions
C This is not in manual? Probably also needs testing

      SUBROUTINE mate_exchange(natoms,off,qoff,pair,qp)
      implicit none

C Inherited integer variables
      integer natoms
C Internal integer variables
      integer i
      integer ca
C Inherited real variables
      real*8 off(3*natoms)
      real*8 qoff(natoms)
      real*8 pair(6*natoms)
      real*8 qp(2*natoms)

c     copy all cations from first string (ca-1 copied!)

      ca = 1

      do i=1,natoms
C For all atoms
         if(qp(i) .gt. 0.d0) then
Copy all but one cation
            qoff(ca) = qp(i)
            off(ca) = pair(i)
            off(ca+natoms) = pair(i+natoms)
            off(ca+2*natoms) = pair(i+2*natoms)
            ca = ca + 1
         end if
      end do

c     copy all anions from second string

      do i=1,natoms
C For all atoms
         if(qp(i+natoms) .lt. 0.d0) then
Copy over all values to new position
            qoff(ca) = qp(i+natoms)
            off(ca) = pair(i+3*natoms)
            off(ca+natoms) = pair(i+4*natoms)
            off(ca+2*natoms) = pair(i+5*natoms)
            ca = ca + 1
         end if
      end do

      return	! Return
      end

C Subroutine mate_1pt_random_sep
C c     Cuts and pastes two selected strings about a single random point
C Alter Z coordinate randomly if max from cluster B is greater than cluster A

      SUBROUTINE mate_1pt_random_sep(natoms,off,pair,idum)
      implicit none

C Inherited integer variables
      integer natoms
      integer idum
C Internal integer variables
      integer i
      integer posit
      integer nvar
C Inherited real variables
      real*8 off(3*natoms)
      real*8 pair(6*natoms)
C Internal real variables 
      real*8 ran3
      real*8 minz
      real*8 maxz
C Import random number generator
      external ran3

      nvar = 3*natoms				! Total coordinates
      posit = ran3(idum)*(natoms-1) + 1		! Choose random atom

C Presumably these commands are run after z_heapsort?
      minz = pair(2*natoms+posit)		! Set to minimum of Cluster A
      maxz = pair(5*natoms+posit+1)		! Set to minimum of cluster B

C Set to zero if cluster A is less than cluster B
      if ( maxz .lt. minz ) minz = 0.d0
      minz = dabs(minz)

C Cut and paste clusters together. Cluster A as part 1...
      do i=1,posit
         off(i) = pair(i)
         off(natoms+i) = pair(natoms+i)
         off(2*natoms+i) = pair(2*natoms+i)
      end do

C And cluster B as part 2...
      do i=posit+1,natoms
         off(i) = pair(i+nvar)
         off(natoms+i) = pair(natoms+i+nvar)
         off(2*natoms+i) = pair(2*natoms+i+nvar) - minz
      end do

      return	! Return
      end
