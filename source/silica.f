c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C Calculator silica structure

      subroutine silica(natoms,na,nb,c1_sisi,c2_sisi,c3_sisi,c4_sisi,
     &                  c1_sio,c2_sio,c3_sio,c4_sio,
     &                  c1_oo,c2_oo,c3_oo,c4_oo,
     &                  ru_sisi,ru_oo,ru_sio,
     &                  xyz,v,diffv,qs,qa,qb,ok)

      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none
      
      real*8 scal
C Integer values received by program
      integer natoms	! Number of atoms
      integer na	! Number of type A
      integer nb	! Number of type B
C Integer values used by program
      integer i
      integer j
C Real values received by program
      real*8 c1_sisi
      real*8 c2_sisi
      real*8 c3_sisi
      real*8 c4_sisi
      real*8 c1_sio
      real*8 c2_sio
      real*8 c3_sio
      real*8 c4_sio
      real*8 c1_oo
      real*8 c2_oo
      real*8 c3_oo
      real*8 c4_oo
      real*8 ru_sisi
      real*8 ru_oo
      real*8 ru_sio
      real*8 xyz(3*natoms)
      real*8 v
      real*8 diffv(3*natoms)
      real*8 qs(natoms)
      real*8 qa
      real*8 qb
C Real values used by program
      real*8 vij
      real*8 dxij
      real*8 dyij
      real*8 dzij
      real*8 rij
      real*8 dv2dr
      real*8 dv2dx(nmax)
      real*8 dv2dy(nmax)
      real*8 dv2dz(nmax)
      real*8 v2tot
      real*8 rij1
      real*8 rhoij1

C Error check
      logical ok

      ok = .true.
      v2tot = 0.d0	! Energy total

C (Original Comments)
c     initialise derivatives

      natoms=na+nb	! Number of atoms

      do i=1,natoms
C For all atoms set derivatives to zero
         dv2dx(i)=0.d0
         dv2dy(i)=0.d0
         dv2dz(i)=0.d0
      end do
      
      do i=1,natoms 
         do j=i+1,natoms
C For all possible combination of atoms

C Work out distance between the two atoms i and j
            dxij=xyz(i)-xyz(j)				! Distance in x
            dyij=xyz(natoms+i)-xyz(natoms+j)		! Distance in y
            dzij=xyz(2*natoms+i)-xyz(2*natoms+j)	! Distance in z
            rij=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)	! Total Distance 

            if (qs(i).eq.qs(j)) then
C If both same atom type...
              if(qs(i).eq.qa) then
C If type is A we have silicon
C Check interaction distance is not to small or to big.
C Probably be nicer if the upper bound was softcoded here
C If there is a problem ok = false, regenerates population
                if (rij.le.ru_sisi.or.rij.gt.20.d0) then
                  ok = .false.
                  return 
                endif

C     (Original Comment)
c     Si-Si interaction
                vij=c1_sisi/rij + c2_sisi/(rij**6) +
     &          c3_sisi*dexp((c4_sisi-rij)/c3_sisi)	! Energy value
                dv2dr=-c1_sisi/(rij**2)-6*c2_sisi/(rij**7)
     &            -dexp((c4_sisi-rij)/c3_sisi)		! Derivative
              else

c     O-O interaction
C Check if unphysical structure
                if (rij.le.ru_oo.or.rij.gt.20.d0) then
                  ok=.false.
                  return
                endif
                vij=c1_oo/rij + c2_oo/(rij**6) +
     &          c3_oo*dexp((c4_oo-rij)/c3_oo)		! Energy value
                dv2dr=-c1_oo/(rij**2)-6*c2_oo/(rij**7)
     &            -dexp((c4_oo-rij)/c3_oo)		! Derivative
              endif
            else
 
c     Si-O interaction
C Check if unphysical structure
              if (rij.le.ru_sio.or.rij.gt.20.d0) then
                ok=.false.
                return
              endif
              vij=c1_sio/rij + c2_sio/(rij**6) +
     &          c3_sio*dexp((c4_sio-rij)/c3_sio)	! Energy Value
              dv2dr=-c1_sio/(rij**2)-6*c2_sio/(rij**7)
     &            -dexp((c4_sio-rij)/c3_sio)		! Derivative
            endif 
     
            v2tot=v2tot+vij		! Add energy to total

C (Original Comment)
c     calculate 2-body first derivatives

C This scale factor appears to be useless
C Maybe it should be softcoded also.
C But currently it achieves nothing
            scal=1.d0     

C Sum up derivatives
            dv2dx(i) = dv2dx(i)+scal*2*dv2dr*dxij/rij	! Add to i...
            dv2dy(i) = dv2dy(i)+scal*2*dv2dr*dyij/rij	!
            dv2dz(i) = dv2dz(i)+scal*2*dv2dr*dzij/rij	!
            dv2dx(j) = dv2dx(j)-scal*2*dv2dr*dxij/rij	! Subtract from j...
            dv2dy(j) = dv2dy(j)-scal*2*dv2dr*dyij/rij	!
            dv2dz(j) = dv2dz(j)-scal*2*dv2dr*dzij/rij	!
            
         end do
      end do
      
      v = v2tot		! Energy is assigned

C (Original comments)
c     calculate total first derivatives
c     diffv is a single 1-d array which stores all derivatives
c     diffv(i)=dv/dx(i), diffv(natoms+i)=dv/dy(i), diffv(2natoms+i)=dv/dz(i)

      do i=1,natoms
C For all atoms, store derivatives and return
         diffv(i)=dv2dx(i)
         diffv(natoms+i)=dv2dy(i)
         diffv(2*natoms+i)=dv2dz(i)
      end do

      return
      end 
