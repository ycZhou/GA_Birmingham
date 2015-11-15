c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C Ionic potential calculator

      subroutine ionic(natoms,B_ab,rho_ab,B_bb,rho_bb,xyz,v,diffv,
     &            qs,qa,qb)

      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none
      
C Inherited integer variables
      integer natoms
C Internal integer values
      integer i
      integer j
C Inherited real variables
      real*8 B_ab
      real*8 rho_ab
      real*8 B_bb
      real*8 rho_bb
      real*8 xyz(3*natoms)
      real*8 v
      real*8 diffv(3*natoms)
      real*8 qs(natoms)
      real*8 qa
      real*8 qb
C Internal Real values
      real*8 velij
      real*8 vbmij
      real*8 dxij
      real*8 dyij
      real*8 dzij
      real*8 rij
      real*8 Bij
      real*8 qij
      real*8 rhoij
      real*8 dv2dr
      real*8 dv2dx(nmax)
      real*8 dv2dy(nmax)
      real*8 dv2dz(nmax)
      real*8 v2tot
      real*8 rij1		! 1/rij
      real*8 rhoij1		! 1/rhoij
 
C     Total energy 
      v2tot = 0.d0

C (Original Comments)
c     initialise derivatives

      do i=1,natoms
C For all atoms set values to zero
         dv2dx(i)=0.d0
         dv2dy(i)=0.d0
         dv2dz(i)=0.d0
      end do
      
      do i=1,natoms
         do j=i+1,natoms
C For all combinations of atoms ...

C Work out distance between i and j
            dxij=xyz(i)-xyz(j)				! Distance apart on x
            dyij=xyz(natoms+i)-xyz(natoms+j)		! Distance apart on y
            dzij=xyz(2*natoms+i)-xyz(2*natoms+j)	! Distance apart on z
            rij=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)	! Total distance (Trig)

C Calculate charge of interaction
            qij=qs(i)*qs(j)

            if (qij .lt. 0.0d0) then
C (Original Comment)
c     A-B interaction
C As we'll have +ve * -ve = -ve
               Bij = B_ab  	!     
               rhoij = rho_ab	!
            else 
               if (qs(i) .lt. 0.d0) then
C If charge on atom A is less than 0
c     B-B interaction
                  Bij = B_bb    	!
                  rhoij = rho_bb	!
               else
c     No A-A Born-Meyer term
                  Bij = 0.0d0		!     
                  rhoij = 1.0d0		!
               end if
            end if
            
C Create reciprocals
            rij1 = 1/rij
            rhoij1 = 1/rhoij

C (Original Comment)
c     electrostatic term
            velij=14.4d0*qij*rij1

C (Original Comment)
c     Born-Meyer term
            vbmij=Bij*dexp(-rij*rhoij1)

C (Original Comment)
c     total PE
            v2tot=v2tot+velij+vbmij	! Current value + electrostatic + Born-Meyer

C (Original Comment)
c     calculate 2-body first derivatives
            dv2dr=(-14.4d0*qij*rij1*rij1)-(Bij*rhoij1)*
     &           dexp(-rij*rhoij1)
            
            dv2dx(i) = dv2dx(i)+dv2dr*dxij*rij1		! Add derivatives to i...
            dv2dy(i) = dv2dy(i)+dv2dr*dyij*rij1		!
            dv2dz(i) = dv2dz(i)+dv2dr*dzij*rij1		!
            dv2dx(j) = dv2dx(j)-dv2dr*dxij*rij1		! And subtract from j...
            dv2dy(j) = dv2dy(j)-dv2dr*dyij*rij1		!
            dv2dz(j) = dv2dz(j)-dv2dr*dzij*rij1		!
            
         end do
      end do
      
C Assign energy to v
      v = v2tot

c     calculate total first derivatives
c     diffv is a single 1-d array which stores all derivatives
c     diffv(i)=dv/dx(i), diffv(natoms+i)=dv/dy(i), diffv(2natoms+i)=dv/dz(i)

      do i=1,natoms
C For all atoms, move derivatives into 1D array
         diffv(i)=dv2dx(i)
         diffv(natoms+i)=dv2dy(i)
         diffv(2*natoms+i)=dv2dz(i)
      end do

      return
      end 
