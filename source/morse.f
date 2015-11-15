c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C Calculate Morse Potential
      SUBROUTINE morse(natoms,alpha,xyz,v,diffv)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Integer value
      integer natoms	! Number of atoms
      integer i		! Counter 1
      integer j		! Counter 2
C Real values
      real*8 alpha		! Potential parameter (Gaussian exponent)
      real*8 xyz(3*natoms) 	! XYZ Coordinates
      real*8 v			! Energy value
      real*8 diffv(3*natoms)
      real*8 aa			! Double alpha
      real*8 afac		! Morse potential (Gaussian function)
      real*8 rij		! Distance between atoms
      real*8 rij1		! Reciprocal of distance
      real*8 dv2dx(nmax)	! Derivatives of x
      real*8 dv2dy(nmax)	! Derivatives of y
      real*8 dv2dz(nmax)	! Derivatives of z
      real*8 dxij		! Difference in x
      real*8 dyij		! Difference in y
      real*8 dzij		! Difference in z
      real*8 dv2dr		! Derivative form on morse potential
C      real*8 dx		! Unused?
C      real*8 dy		! Unused?
C      real*8 dz		! Unused?
      real*8 dvdx		! First derivative for current atom in x
      real*8 dvdy		! First derivative for current atom in y
      real*8 dvdz		! First derivative for current atom in z

C     (Original Comments)
c     Calculates the energy and first derivatives of the energy
c     of the cluster in the xyz array
c     
              
      aa=2.0d0*alpha
      v = 0.d0

c ... initialise derivatives
      
      do i=1,natoms
         dv2dx(i)=0.d0
         dv2dy(i)=0.d0
         dv2dz(i)=0.d0
      end do
      
C For all atoms except the last
      do i=1,natoms-1
C Compare to all atoms except the first         
         do j=i+1,natoms
            
C Calculate atom displacements between each other
C This is just taken from simple trigonometry

            dxij=xyz(i)-xyz(j)			! Difference in x coordinates
            dyij=xyz(natoms+i)-xyz(natoms+j)	! ... in y coordinates
            dzij=xyz(2*natoms+i)-xyz(2*natoms+j)	! ... in z coordinates
            rij=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)	! Distance between atoms
            
C Calculate exponential Gaussian term
C This is of the form e^(alpha(1-[distance between i and j]))

	    afac=dexp(alpha*(1.0d0-rij))

C Energy value
C This appears to be a summation of the above Gaussian terms, close to squared 
C in the form afac*(afac-2) i.e. 2 away from being afac^2
            v=v+afac*(afac-2.0d0)

C (Original comment)
c ... calculate 2-body first derivatives

C Value of derivative of morse potential term
C Of the form 2(alpha)(gaussian term)(1-gaussian term)
C This can easily be reproduced analytically if we want

            dv2dr=aa*afac*(1.0d0-afac)

C Reciprocal of distance between r and i
            rij1 = 1.d0/rij

C From this we can work out the values for x, y and z
C as shown here they are derivative of equation * distance in plane * reciprocal distance
C For the current pair of atoms

            dvdx = dv2dr*dxij*rij1
            dvdy = dv2dr*dyij*rij1
            dvdz = dv2dr*dzij*rij1

C And this can be summed for the current atom as positive forces on
C atom a

            dv2dx(i) = dv2dx(i) + dvdx
            dv2dy(i) = dv2dy(i) + dvdy
            dv2dz(i) = dv2dz(i) + dvdz

C And negative forces on atom b

            dv2dx(j) = dv2dx(j) - dvdx
            dv2dy(j) = dv2dy(j) - dvdy
            dv2dz(j) = dv2dz(j) - dvdz

         end do
         
      end do

C (Original Comment)
c ... total potential energy

c ... calculate total first derivatives
c ... diffv is a single 1-d array which stores all derivatives
c ... diffv(i)=dv/dx(i), diffv(natoms+i)=dv/dy(i), diffv(2*natoms+i)=dv/dz(i)
C No further explanation required

      do i=1,natoms
         diffv(i)=dv2dx(i)
         diffv(natoms+i)=dv2dy(i)
         diffv(2*natoms+i)=dv2dz(i)
      end do
      
      return
      end
