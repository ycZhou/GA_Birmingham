c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C This is the Murrell Mottram potential calculator
      SUBROUTINE murrmott(natoms,idamp,ihard,de,re,a2,a3,coeff,
     &            rcut2,xyz,v,diffv)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none

C Firstly lets look at incoming variables
C Integer values for Murrell Motram
      integer natoms	! Number of atoms
      integer idamp	! Dampening function (1 = exp, 2 = tanh, 3 = gaussian, 4 = sech)
      integer ihard	!
C Real values
      real*8 de		! Scale to eV
      real*8 re
      real*8 a2
      real*8 a3
      real*8 coeff(11)		! Coefficients
      real*8 rcut2		! Cutoff distance squared
      real*8 xyz(3*natoms)	! xyz coordinates
      real*8 v			! Energy
      real*8 diffv(3*natoms)	!
      
C Now we have internal (dynamic) variables
C Integer values
      integer i		! Atom a
      integer j		! Atom b
      integer k		! Atom c

C Real values
      real*8 dv2dx(nmax)	! Derivatives with respect to x
      real*8 dv2dy(nmax)	! Derivatives with respect to y
      real*8 dv2dz(nmax)	! Derivatives with respect to z
      real*8 tr
      real*8 a22		! a2^2
      real*8 dv3dx(nmax)	! Derivatives of i-k with respect to x
      real*8 dv3dy(nmax)	! Derivatives of i-k with respect to y
      real*8 dv3dz(nmax)	! Derivatives of i-k with respect to z
      real*8 v2tot		! Total Energy with respect to i-j
      real*8 v3tot		! Total Energy with respect to i-k
      real*8 dxij		! Distance for x between i and j
      real*8 dyij		! Distance for y between i and j
      real*8 dzij		! Distance for z between i and j
      real*8 rij		! Distance between atoms 1 and 2
      real*8 rhoij		! rij - 1
      real*8 a2r		! a2*(rij - 1)
      real*8 exar		! Contributions to derivative of potential
      real*8 ex50		! Contributions to derivative of potential
      real*8 dxik		! Distance for x between i and k
      real*8 dyik		! Distance for y between i and k
      real*8 dzik		! Distance for z between i and k
      real*8 dxjk		! Distance for x between j and k
      real*8 dyjk		! Distance for y between j and k
      real*8 dzjk		! Distance for z between j and k
      real*8 rik		! Distance between atoms 1 and 3 (i and k)
      real*8 rjk		! Distance between atoms 2 and 3
      real*8 rhoik		! rik - 1 
      real*8 rhojk		! rjk - 1
      real*8 q1			! Q coordinate (x)
      real*8 q2			! Q coorindate (y)
      real*8 q3			! Q coordinate (z)
      real*8 q11		! Q(x)^2
      real*8 q22		! Q(y)^2
      real*8 q33		! Q(z)^2
      real*8 damp
      real*8 ddamp
      real*8 t
      real*8 p
      real*8 dpdq1
      real*8 dpdq2
      real*8 dpdq3
      real*8 v3
      real*8 dv3dq1
      real*8 dv3dq2
      real*8 dv3dq3
      real*8 dv3rij
      real*8 dv3rik
      real*8 dv3rjk
      real*8 dvdxij
      real*8 dvdxik
      real*8 dvdxjk
      real*8 dvdyij
      real*8 dvdyik
      real*8 dvdyjk
      real*8 dvdzij
      real*8 dvdzik
      real*8 dvdzjk
      real*8 fhard
      real*8 dv2dr
      real*8 rij2		! (Distance between 1 and 2)^2
      real*8 rik2               ! (Distance between 1 and 3)^2
      real*8 rjk2               ! (Distance between 2 and 3)^2
      real*8 rij1		! Reciprocal of rij
      real*8 rik1		! Reciprocal of rik
      real*8 rjk1		! Reciprocal of rjk
      real*8 dvdx
      real*8 dvdy
      real*8 dvdz
      real*8 g1			! Constants
      real*8 g2			!
      real*8 g3			!		
      real*8 g4			!
      real*8 g5			!
      real*8 g6			!
      real*8 g7			!
      real*8 g8			!
      real*8 c0			!
      real*8 c1			!
      real*8 c2			!
      real*8 c3			!
      real*8 c4			!
      real*8 c5			!
      real*8 c6			!
      real*8 c7			!
      real*8 c8			!
      real*8 c9			!
      real*8 c10		!
C Boolean variable
      logical flag
C Static variables
      save flag
C Irreducible representations of S3 permutation group
      save g1		! sqrt(1/3)
      save g2		! sqrt(1/2)
      save g3		! sqrt(1/6)
      save g4		! 1/3
      save g5		! 1/6
      save g6		! sqrt(2/3)
      save g7		! sqrt(2)/3
      save g8		! sqrt(1/18)
      save tr		!
      save a22		! a2 squared
      save fhard	!
      save c0		! Coefficients 1...
      save c1		!..
      save c2		!.. 
      save c3		!..
      save c4		!..
      save c5		!..
      save c6		!..
      save c7		!..
      save c8		!..
      save c9		!..
      save c10		!.. to 11
C Set flag value
      data flag /.true./

      if (flag) then

C (Original Comments)
c ... define constants

         g1=dsqrt(1.0d0/3.0d0)
         g2=dsqrt(0.5d0)
         g3=dsqrt(1.0d0/6.0d0)
         g4=1.d0/3.d0
         g5=1.d0/6.d0
         g6=dsqrt(2.d0/3.d0)
         g7=dsqrt(2.d0)/3.d0
         g8=dsqrt(1.d0/18.d0)
         tr=(1.d0+dsqrt(5.d0))/2.d0
         a22=a2*a2
         fhard=dble(ihard)

C Extract Coefficients

         c0 = coeff(1)
         c1 = coeff(2)
         c2 = coeff(3)
         c3 = coeff(4)
         c4 = coeff(5)
         c5 = coeff(6)
         c6 = coeff(7)
         c7 = coeff(8)
         c8 = coeff(9)
         c9 = coeff(10)
         c10 = coeff(11)
C Static variables defined, therefore make flag false
         flag = .false.

      end if

C Zero energy totals
   
      v2tot = 0.d0
      v3tot = 0.d0

C (Original Comments)
c ... initialise derivatives

      do i=1,natoms
         dv2dx(i)=0.d0
         dv2dy(i)=0.d0
         dv2dz(i)=0.d0
         dv3dx(i)=0.d0
         dv3dy(i)=0.d0
         dv3dz(i)=0.d0
      end do

C For all atoms except the last
      do i=1,natoms-1
C Compare to all atoms except the current
         do j=i+1,natoms
      
            dxij=xyz(i)-xyz(j)			! Calculate distance difference in x-axis
            dyij=xyz(natoms+i)-xyz(natoms+j)	! Calculate distance difference in y-axis
            dzij=xyz(2*natoms+i)-xyz(2*natoms+j)	! Calculate distance difference in z-axis
            rij2=dxij*dxij+dyij*dyij+dzij*dzij	! Calculate total distance, squared

C If this squared value is less than the cutoff squared
            if (rij2 .le. rcut2) then
C Calculate distance between i and j
               rij=dsqrt(rij2)
C This gives our distance in reduced coordinates
C Rho is the (distance between i and j)-1
C               rhoij=rij-1.d0
C And this gives rho is real coordinates
               rhoij=(rij/re)-1

               a2r=a2*rhoij	! And thus a2r is this rho times a2

C (Original Comments)
c       ... calculate 2-body potential (including hard wall)
C So here we are calculating the energy value. Bit of a mouthful...                 

               exar=dexp(-a2r)				! = e^(-a2*(rij-1))
               ex50=fhard*dexp(-50.d0*(rhoij+0.2d0))	! = fhard*e^(-50*(rij-0.8))
               v2tot=v2tot-(1.d0+a2r)*exar+ex50 	! running total = total-[e^(-a2*(rij-1))]*(1+[a2*(rij-1)]+[fhard*e^(-50*(rij-0.8))]
      
C (Original Comments)
c       ... calculate 2-body first derivatives
C And here we calculate the first derivatives                 

C              First Derivative
               dv2dr=a22*rhoij*exar-50.d0*ex50		! a2^2 * (rij-1) * e^(-a2*(rij-1)) - ( 50 * fhard*e^(-50*(rij-0.8)))

               rij1 = 1.d0/rij				! reciprocal of rij

               dvdx = dv2dr*dxij*rij1			! Derivative of x component
               dvdy = dv2dr*dyij*rij1			! Derivative of y component
               dvdz = dv2dr*dzij*rij1			! Derivative of z component
      
C              Add to sum of derivatives for i
               dv2dx(i)=dv2dx(i)+dvdx
               dv2dy(i)=dv2dy(i)+dvdy 
               dv2dz(i)=dv2dz(i)+dvdz
C              Subtract from sum of derivatives for j
               dv2dx(j)=dv2dx(j)-dvdx
               dv2dy(j)=dv2dy(j)-dvdy
               dv2dz(j)=dv2dz(j)-dvdz

C       Now to compute the 3 body interaction
C       (Original Comment)
c       ... 3-body loop
      
               do k=j+1,natoms

C (Original comment)
c ... set up triangle
 
C       Now configure interactions between atoms 1 and 3
                  dxik=xyz(i)-xyz(k)				! Distance between x coordinates for atom i and k
                  dyik=xyz(natoms+i)-xyz(natoms+k)		! Distance between y coordinates for atom i and k
                  dzik=xyz(2*natoms+i)-xyz(2*natoms+k)		! Distance between z coordinates for atom i and k
C       Configure interactions between atoms 2 and 3
                  dxjk=xyz(j)-xyz(k)				! Distance between x coordinates for atom j and k
                  dyjk=xyz(natoms+j)-xyz(natoms+k)		! Distance between y coordinates for atom j and k
                  dzjk=xyz(2*natoms+j)-xyz(2*natoms+k)		! Distance between z coordinates for atom j and k
                  rik2=dxik*dxik+dyik*dyik+dzik*dzik		! Distance (i to k)^2
                  rjk2=dxjk*dxjk+dyjk*dyjk+dzjk*dzjk		! Distance (j to k)^2

C If distances squared are less than cutoff for 3-body interaction...                  
                  if (rik2 .le. rcut2 .or. rjk2 .le. rcut2) then
                     
                     rik=dsqrt(rik2)				! Distance between i and k
                     rjk=dsqrt(rjk2)				! Distance between j and k

C                    This will give our result in reduced coordinates                     
C                     rhoik=rik-1.d0				! rik -1
C                     rhojk=rjk-1.d0				! rjk -1

C                    And this gives our values in real coordinates 
                     rhoik=(rik/re)-1                            ! rik -1
                     rhojk=(rjk/re)-1                            ! rjk -1
                    
                     q1=g1*(rhojk+rhoij+rhoik)			! g1*rij + g1*rjk + g1*rik - 3*g1 
                     q2=g2*(rhoij-rhoik)			! g2*rij - g2*rik + 2*g2
                     q3=g3*(2.d0*rhojk-rhoij-rhoik)		! 3*g3 - g3*rij - 2*g3*rjk - g3*rik
                     
                     q11=q1*q1					! Q(x)^2
                     q22=q2*q2					! Q(y)^2
                     q33=q3*q3					! Q(z)^2
                     
C (Original Comments)
c ... calculate 3-body potential v3 = p(q1,q2,q3)damp(q1)
c ... decaying part damp(q1) and its 1st derivs ddamp(q1)

                     if (idamp.eq.1)then
c     ... exponential
                        damp=dexp(-a3*q1)
                        ddamp=-a3*damp
                        
                     else if (idamp.eq.2) then
c     ... tanh
                        t=tanh(a3*q1*0.5d0)
                        damp=(1.0d0-t)*0.5d0
                        ddamp=a3*0.25d0*(t*t-1.0d0)
                        
                     else if (idamp.eq.3) then
c     ... gaussian
                        damp=dexp(-a3*q1*q1)
                        ddamp=-2.d0*a3*q1*damp
                        
C     ... sech
                     else if (idamp.eq.4) then

                        damp=(1.0d0/cosh(a3*q1))
                        t=tanh(a3*q1)
                        ddamp=-a3*t*damp
                        
                     end if

C (Original Comments)
c ... polynomial p(q1,q2,q3)

                     p=c0+q1*(c1+c2*q1+c4*q11+c7*q1*q11)
     &                    +(q22+q33)*(c3+c5*q1+c8*q11+c9*(q22+q33))
     &                    +q3*(q33-3.d0*q22)*(c6+c10*q1)

C (Original Comments)
c ... 1st derivs of p with respect to Q
                       
                     dpdq1=c1+q1*(2.d0*c2+3.d0*c4*q1+4.d0*c7*q11)
     &                    +(q22+q33)*(c5+2.d0*c8*q1)
     &                    +c10*q3*(q33-3.d0*q22)
     
                     dpdq2=2.d0*q2*(c3+c5*q1+c8*q11+2.d0*c9*(q22+q33)
     &                    -3.d0*q3*(c6+c10*q1))
                       
                     dpdq3=2.d0*q3*(c3+c5*q1+c8*q11+2.d0*c9*(q22+q33))
     &                    +3.d0*(q33-q22)*(c6+c10*q1)

C (Original comments)
c       ... 3-body potential v3

                     v3=p*damp			! Polynomial * Dampening
                     v3tot=v3tot+v3		! Total Energy

c       ... 1st derivs of v3
                    
                     dv3dq1=damp*dpdq1+p*ddamp	! Q1
                     dv3dq2=damp*dpdq2		! Q2
                     dv3dq3=damp*dpdq3		! Q3

 
                     dv3rij=g1*dv3dq1+g2*dv3dq2-g3*dv3dq3	! 3-body derivative with respect to ij
                     dv3rik=g1*dv3dq1-g2*dv3dq2-g3*dv3dq3	! 3-body derivative with respect to jk (Same as above)
                     dv3rjk=g1*dv3dq1+2.d0*g3*dv3dq3		! 3-body derivative with respect to jk

                     rik1 = 1.d0/rik		! Reciprocal for derivative
                     rjk1 = 1.d0/rjk		! Reciprocal for derivative
                     
C 		     Calculate derivatives of x, y and z for ij...
                     dvdxij=dv3rij*dxij*rij1
                     dvdyij=dv3rij*dyij*rij1
                     dvdzij=dv3rij*dzij*rij1
C 	             Calculate derivatives of x, y and z for ik 
                     dvdxik=dv3rik*dxik*rik1
                     dvdyik=dv3rik*dyik*rik1
                     dvdzik=dv3rik*dzik*rik1
C 		     And calculate derivatives of x, y and z for jk
                     dvdxjk=dv3rjk*dxjk*rjk1
                     dvdyjk=dv3rjk*dyjk*rjk1
                     dvdzjk=dv3rjk*dzjk*rjk1

C                    Sum derivatives for each atom, starting with i...
                     dv3dx(i)=dv3dx(i)+dvdxij+dvdxik
                     dv3dy(i)=dv3dy(i)+dvdyij+dvdyik
                     dv3dz(i)=dv3dz(i)+dvdzij+dvdzik
C                    Then subtracting from sum of j (addinging 3rd body influence) 
                     dv3dx(j)=dv3dx(j)-dvdxij+dvdxjk
                     dv3dy(j)=dv3dy(j)-dvdyij+dvdyjk
                     dv3dz(j)=dv3dz(j)-dvdzij+dvdzjk
C                    And then subtracting derivatives for all parts of k
                     dv3dx(k)=dv3dx(k)-dvdxik-dvdxjk
                     dv3dy(k)=dv3dy(k)-dvdyik-dvdyjk
                     dv3dz(k)=dv3dz(k)-dvdzik-dvdzjk
                     
                  end if
                    
               end do

            end if
            
         end do
         
      end do
      
C (Original Comment)
c     ... scale energies (units = eV)
      
      v2tot=v2tot*de
      v3tot=v3tot*de
     
C (Original Comment)
c     ... total potential energy
c     NB during optimisation procedure scale energy and derivatives by 
c     fscale to speed up convergence
c     and avoid memory problems for large clusters

      v=v2tot+v3tot	! Sum energies

C (Original Comment)
c     ... calculate total first derivatives
c     ... diffv is a single 1-d array which stores all derivatives
c     ... diffv(i)=dv/dx(i), diffv(nat+i)=dv/dy(i), diffv(2nat+i)=dv/dz(i)

      do i=1,natoms
         diffv(i)=(dv2dx(i)+dv3dx(i))*de		! Combine derivatives for x and scale
         diffv(natoms+i)=(dv2dy(i)+dv3dy(i))*de		! Combine derivatives for y and scale
         diffv(2*natoms+i)=(dv2dz(i)+dv3dz(i))*de	! Combine derivatives for z and scale
      end do

      return
      end
