c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C Generate population for silica
      subroutine gen_pop_silica(na,nb,nclust,nstride,idum,pop,q,qa,qb,
     &           ru_sisi,ru_sio,ru_oo,range)
      implicit none

C Integer Values
      integer na
      integer nb
      integer nclust
      integer nstride
      integer idum
      integer i
      integer j
C Real Values
      real*8 pop((3*(na+nb)+1)*nclust)
      real*8 q((na+nb)*nclust)
      real*8 qa
      real*8 qb
      real*8 ru_sisi
      real*8 ru_sio
      real*8 ru_oo
      real*8 range
      real*8 ran3
C Boolean check
      logical ok
C Import random number generator
      external ran3

C (Original comments)
c     seed random number generator and generate initial population
       
C Start point
       i=1
C Create loop which runs until complete
      do 
         if (i.le.nclust) then
C Create coordinates for atoms in cluster
         do j=1,3*(na+nb)
            pop(j+(i-1)*nstride) = 3*ran3(idum)*range
         end do
C Create charges for atoms in cluster
         do j=1,na
            q(j+(i-1)*(na+nb)) = qa
         end do
         do j=na+1,na+nb
            q(j+(i-1)*(na+nb)) = qb
         end do
c     reject unphysical clusters using boolean flag
         ok=.true.
C Check if cluster is physically possible
         call check_unphys_silica(ok,na,nb,pop(1+(i-1)*nstride),
     &           q(1+(i-1)*(na+nb)),qa,qb,ru_sisi,ru_sio,ru_oo)
         if (ok) i=i+1
         end if
C Exit if i is greater than number of clusters
         if (i.eq.nclust+1) exit
      end do

      return
      end

C Subroutine to check physicality of Si clusters
      subroutine check_unphys_silica(ok,na,nb,xyz,qs,qa,qb,
     &           ru_sisi,ru_sio,ru_oo)
      implicit none

C Integer values
      integer na	! Number of a
      integer nb	! Number of b
      integer i
      integer j
C Real Values
      real*8 xyz(3*(na+nb))
      real*8 qs(na+nb)
      real*8 qa		! Charge of a
      real*8 qb		! Charge of b
      real*8 ru_sisi
      real*8 ru_sio
      real*8 ru_oo
      real*8 range	! Atom range
      real*8 dxij 	! X Distance between i and j
      real*8 dyij	! Y Distance between i and j
      real*8 dzij	! Z Distance between i and j
      real*8 rij	! Radial distance between i and j
C Boolean Flags
      logical ok

C Loop over all atoms
      do i=1,na+nb  
C Loop over all other atoms
         do j=i+1,na+nb  
C Calculate distance between the two
            dxij=xyz(i)-xyz(j)
            dyij=xyz(na+nb+i)-xyz(na+nb+j)
            dzij=xyz(2*(na+nb)+i)-xyz(2*(na+nb)+j)
C Some simple trigonometry for you :)
C Now check if this is within our parameters
            rij=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)
            if (qs(i).eq.qs(j)) then
              if(qs(i).eq.qa) then
                 if (rij.le.ru_sisi) ok=.false. 
              else
                 if (rij.le.ru_oo) ok=.false.
              endif
            else
              if (rij.le.ru_sio) ok=.false.
            endif
         enddo
      enddo 
      return
      end
C And we are done
