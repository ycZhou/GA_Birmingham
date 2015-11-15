c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C This is the Gupta potential calculator

      subroutine gupta(natoms,nspec,xyz,qoff,f,g,aij,vij,pij,qij
     &     ,r0,iflag,cutoff_gupta,cutoff_start,cutoff_end,a3a,a4a,a5a,
     &     x3x,x4x,x5x)

C     Import variables
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none
      
C     Integer values for ga
      integer natoms	! Number of atoms
      integer nspec	! Number of species
      integer i		! Atom A reference
      integer j		! Atom B reference
      integer ii	! Access to parameters array
      integer nn	! Counter for array access of derivatives

C     Real values for calculation imported from main program
      real*8 xyz(3*natoms)	! Coordinate systems
      real*8 qoff(natoms)	! Charge on systems
      real*8 f   		!
      real*8 g(3*natoms)	!
      real*8 aij(nspec*nspec)	! Gupta parameters: A
      real*8 vij(nspec*nspec)	! zeta
      real*8 pij(nspec*nspec)	! p
      real*8 qij(nspec*nspec)	! q
      real*8 r0(nspec*nspec)	! bond length
      real*8 cutoff_start(nspec*nspec)	! Coefficients from file: cutoff start
      real*8 cutoff_end(nspec*nspec)	! cutoff end
      real*8 a3a(nspec*nspec)	! a3
      real*8 a4a(nspec*nspec)	! a4
      real*8 a5a(nspec*nspec)	! a5
      real*8 x3x(nspec*nspec)	! x3
      real*8 x4x(nspec*nspec)	! x4
      real*8 x5x(nspec*nspec)	! x5

C     Real values used internally
      real*8 e2			!
      real*8 rhoij		! Gupta parameters for current variable: p
      real*8 aaij		! A
      real*8 z2ij		! zeta
      real*8 r0ij		! Bond length
      real*8 rdiv		! (rij * r0) -1
      real*8 qqij		! q
      real*8 v			! Total Energy
      real*8 v1(nmax)		! Repulsive Arrays
      real*8 v2(nmax)		! Attractive Arrays
      real*8 v1ij		! Repulsive Term
      real*8 v2ij		! Attractive Term
      real*8 v1i		! Repulsive Value for Atom A
      real*8 v2i		! Attractive Value for Atom B
      ! nspecmax2=9
      real*8 c1ij(9)	! Term used in derivative calculation
      real*8 c2ij(9)	! Term used in derivative calculation
      real*8 rij		! Distance between A and B
      real*8 dxij		! Distance between A and B in x
      real*8 dyij		! Distance between A and B in y
      real*8 dzij		! Distance between A and B in z
      real*8 dv1dx(nmax)	! Derivatives of Repulsive Term (x axis)
      real*8 dv1dy(nmax)	! Derivatives of Repulsive Term (y axis)
      real*8 dv1dz(nmax)	! Derivatives of Repulsive Term (z axis)
      real*8 v2i1               ! Reciprocal of Attractive Term
      real*8 dv2dx(nmax)	! Derivative of Attractive Term (x axis)
      real*8 dv2dy(nmax)	! Derivative of Attractive Term (x axis)
      real*8 dv2dz(nmax)	! Derivative of Attractive Term (x axis)
      real*8 dv2dxij(nmax**2/2)	! Total Array of Attractive derivatives for x
      real*8 dv2dyij(nmax**2/2)	! Total Array of Attractive derivatives for y
      real*8 dv2dzij(nmax**2/2)	! Total Array of Attractive derivatives for z
      real*8 dv1dxi		! Repulsive running total of derivative x for atom A
      real*8 dv1dyi		! Repulsive running total of derivative y for atom A
      real*8 dv1dzi		! Repulsive running total of derivative z for atom A
      real*8 dv2dxi		! Repulsive running total of derivative x for atom B
      real*8 dv2dyi		! Repulsive running total of derivative y for atom B
      real*8 dv2dzi		! Repulsive running total of derivative z for atom B
      real*8 dv1dxij		!
      real*8 dv1dyij		!
      real*8 dv1dzij		!
      real*8 cc1ij		!
      real*8 cc2ij		!
      real*8 rij1		! reciprocal of rij
      real*8 xval		!
      real*8 yval		!
      real*8 zval		!
      real*8 repulsive		! Derivative for A B of repulsive term
      real*8 attractive		! Derivative for A B of attractive term

C ADDED BY ANDREW LOGSDAIL, OCTOBER 2009
C For use with gupta cutoff

      real*8 use_cs	! Cutoff start
      real*8 use_ce	! Cutoff end
      real*8 use_a3a	! Coefficients: a3
      real*8 use_a4a	! .. a4
      real*8 use_a5a	! .. a5
      real*8 use_x3x	! .. x3
      real*8 use_x4x	! .. x4
      real*8 use_x5x	! .. x5
      real*8 factor1	! .. rij - ce
      real*8 factor2	! .. (rij - ce)^2
      real*8 factor3	! .. (rij - ce)^3
      real*8 factor4	! .. (rij - ce)^4
      real*8 factor5	! .. (rij - ce)^5

C     Boolean operator implying cutoff use
      logical cutoff_gupta
C     Other boolean operators
      logical flag		! Set parameters
      logical iflag		! Check for errors

C     Set constant value
      data flag /.true./

C     Make constants
      save flag
      save c1ij
      save c2ij

C     Set cutoff calues
      use_cs=0.d0
      use_ce=0.d0

C     Initiate variables on first run
      if (flag) then
         
C     For all types on atoms species
         do i=1,nspec*nspec
            r0(i) = 1.d0/r0(i)		! Create reciprocal of bond length
            vij(i) = vij(i)*vij(i)	! Square zeta
            c1ij(i) = -pij(i)*r0(i) 	! Negative repulsive value * reciprocal of bond length
            c2ij(i) = qij(i)*r0(i)	! Attractive term * reciprocal of bond length
         end do
         
C     Variables initiated. Do not repeat
         flag = .false.
      end if

C     For all atoms, initiate values on this run
      do i=1,natoms
         v1(i) = 0.d0
         v2(i) = 0.d0
         dv1dx(i) = 0.d0	! Differential of x
         dv1dy(i) = 0.d0	! Differential of y
         dv1dz(i) = 0.d0	! Differential of z
         dv2dx(i) = 0.d0 	! Second Derivative of x
         dv2dy(i) = 0.d0	! Second Derivative of y
         dv2dz(i) = 0.d0	! Second Derivative of z
      end do

      iflag = .false.		! Set flag for completion
      v = 0.d0			! Set energy to zero
      nn = 0			! Counter

C For all atoms...
      do i=1,natoms

C Set values to those previously calculated
         v1i = v1(i)		! Set Repulsive to last value
         v2i = v2(i)		! Set Attractive to last value
         dv1dxi = dv1dx(i)	! Set derivative as value from last x
         dv1dyi = dv1dy(i)	! Set derivative as value from last y
         dv1dzi = dv1dz(i)	! Set derivative as value from last z
         dv2dxi = dv2dx(i)	! Set 2nd derivative as value from last x
         dv2dyi = dv2dy(i)	! Set 2nd derivative as value from last y
         dv2dzi = dv2dz(i)	! Set 2nd derivative as value from last z

C Then let us pair up with all other atom types not previously paired
         do j=i+1,natoms

C (Original Comment)
c     get correct potential parameters
            nn = nn + 1			! Set counter position

C           Check number of species to define access to array
            if (nspec .eq. 2) then
               ii = qoff(j) + (qoff(i) - 1) * nspec
            else
               ii=1
            end if

C           Get values
            r0ij = r0(ii)		! bond length
            rhoij = pij(ii)		! p 
            aaij = aij(ii)		! A
            z2ij = vij(ii)		! zeta
            qqij = qij(ii)		! q
            cc1ij = c1ij(ii)		! Repulsive term for use with derivatives
            cc2ij = c2ij(ii)		! Attractive term for use with derivatives

c     get correct cutoff parameters if used

            if (cutoff_gupta) then
               use_cs=cutoff_start(ii)	! Cs
               use_ce=cutoff_end(ii)	! Ce
               use_a3a=a3a(ii)		! a3
               use_a4a=a4a(ii)		! a4
               use_a5a=a5a(ii)		! a5
               use_x3x=x3x(ii)		! x3
               use_x4x=x4x(ii)		! x4
               use_x5x=x5x(ii)		! x5
            end if

C (Original Comments)
c     calculate exponentials

C Firstly let us calculate the distance between atoms
            dxij = xyz(i) - xyz(j)			! Distance in x axis
            dyij = xyz(i+natoms) - xyz(j+natoms)	! Distance in y axis
            dzij = xyz(i+2*natoms) - xyz(j+2*natoms)	! Distance in z axis
            rij = dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)	! Total Distance

C (Original comments)
c     dodgy "if" - some small clusters seem to explode
c     setting rij small keeps gradients small -> quicker convergence
C Added comment: Would be good to softcode this as may prevent big clusters from
C                completing calculation. Currently value of 100A (10nm)

            if (rij .gt. 100.d0) goto 99	! Cluster has exploded, create new
            rdiv = rij*r0ij - 1.d0		! (rij * r0) - 1
            rij1 = 1.d0/rij			! Reciprocal of rij

C (Original Comments)
c     potential

C     Time to add the if statements for the gupta potential.
C     Andrew Logsdail 2009

            attractive = 0.d0	! Derivative Attractive Term
            repulsive = 0.d0	! Derivative Repulsive Term

            if (.not. cutoff_gupta .or. rij .lt. use_cs) then
C           If distance is less than cutoff, or cutoff not being used

              v1ij = aaij * exp(-rhoij*rdiv)		! Repulsive Term
              v2ij = z2ij * exp(-2.d0*qqij*rdiv)	! Attractive Term

C             Calculate derivative form
              repulsive = v1ij * cc1ij	! Repulsive * -p * 1/r0
              attractive = v2ij * cc2ij	! Attractive * q * 1/r0

            else if (rij .lt. use_ce) then
C           If distance between A and B is less than the end of the cutoff distance

C             Calculate factors for polynomial cutoff
              factor1 = rij-use_ce		!
              factor2 = factor1*factor1		!
              factor3 = factor2*factor1		!
              factor4 = factor3*factor1		!
              factor5 = factor4*factor1		!

C             Define attractive and repulsive components as 5th order polynomials
              v1ij = use_a5a*factor5 + use_a4a*factor4 + use_a3a*factor3
              v2ij = use_x5x*factor5 + use_x4x*factor4 + use_x3x*factor3

C             Derivative of attractive term
              attractive = -((5*use_x5x*factor4 + 4*use_x4x*factor3 +
     & 3*use_x3x*factor2)*v2ij)
C             Derivative of repulsive term
              repulsive = 5*use_a5a*factor4 + 4*use_a4a*factor3 +
     & 3*use_a3a*factor2

C             Square attractive component as this will be square rooted later.
              v2ij = v2ij*v2ij

            else

C             Beyond cutoff set values to zero
              v1ij = 0
              v2ij = 0

            end if

C           Add new values for attraction/repulsion to previous for atom A
            v1i = v1i + v1ij
            v2i = v2i + v2ij

C           Store attraction/repulsion for Atom B in array
            v1(j) = v1(j) + v1ij
            v2(j) = v2(j) + v2ij

C (Original Comment)
c     derivatives of V1

C           Calculate derivatives for each axis
            dv1dxij = repulsive*dxij*rij1	! Repulsive * dxij * 1/rij
            dv1dyij = repulsive*dyij*rij1	! Repulsive * dyij * 1/rij
            dv1dzij = repulsive*dzij*rij1	! Repulsive * dzij * 1/rij

C           Derivatives added to previous values for Atom I
            dv1dxi = dv1dxi + dv1dxij
            dv1dyi = dv1dyi + dv1dyij
            dv1dzi = dv1dzi + dv1dzij

C           And subtracted from previous values for Atom J
            dv1dx(j) = dv1dx(j) - dv1dxij
            dv1dy(j) = dv1dy(j) - dv1dyij
            dv1dz(j) = dv1dz(j) - dv1dzij

C (Original Comments)
c     derivatives of V2
            
c     common quantity
c     qia/r0ia * via * xai/rai
c     c2ia = qia/r0ia 

C Not sure what these values ^^^^^^ are for
C As they are not used in the program. 
C We will presume it is redundant code

C           Calculate derivatives of attractive term for each axis and store
            dv2dxij(nn) = attractive*dxij*rij1 
            dv2dyij(nn) = attractive*dyij*rij1 
            dv2dzij(nn) = attractive*dzij*rij1 

C (Original Comment)
c     accumulate sum qia/r0ia * via * xai/rai

C           Attractive Derivatives added to previous values for Atom I
            dv2dxi = dv2dxi + dv2dxij(nn)
            dv2dyi = dv2dyi + dv2dyij(nn)
            dv2dzi = dv2dzi + dv2dzij(nn)

C           Attractive Derivatives subtracted from previous values for Atom J
            dv2dx(j) = dv2dx(j) - dv2dxij(nn)
            dv2dy(j) = dv2dy(j) - dv2dyij(nn)
            dv2dz(j) = dv2dz(j) - dv2dzij(nn)

         end do
C So we have looped over all atoms j, and calculated all derivatives
C and energy values given the current arrangement of atoms with respect to atom i
C Now let us review the energy values for this single atom

C        Reciprocal of sqrt of attractive total for atom i
         v2i1 = 1.d0/dsqrt(v2i)

C        Total energy is current total plus + repulsive term - attractive term
C        In keeping with the Gupta potential form
         v = v + ( v1i - dsqrt(v2i))

C (Original Comment)
c     this should be ok since v2(i) is not used again in this loop
C ? We are saving the reciprocal of the attractive term for this atoms
         v2(i) = v2i1

C Save derivative of repulsive to array
         dv1dx(i) = dv1dxi
         dv1dy(i) = dv1dyi
         dv1dz(i) = dv1dzi

C Save (derivative of attractive for atom i) * (1/Attractive Energy) to array
         dv2dx(i) = v2i1 * dv2dxi
         dv2dy(i) = v2i1 * dv2dyi
         dv2dz(i) = v2i1 * dv2dzi

      end do
C At this point we have looped interactions between all atom types i and j

C (Original Comment)
c     sum 2nd part of derivatives

C (Reset value which maps across array with all attractive terms
      nn = 0.d0

      do i=1,natoms
         do j=i+1,natoms
C For all atoms, pair with all others and compare...

            nn = nn + 1			! Array pointer

C	    Attractive derivative for j is (current value) - (reciprocal of ij energy) * (derivative)
            dv2dx(j) = dv2dx(j) - v2(i)*dv2dxij(nn)
            dv2dy(j) = dv2dy(j) - v2(i)*dv2dyij(nn)
            dv2dz(j) = dv2dz(j) - v2(i)*dv2dzij(nn)
 
C           Attractive derivative for i is (current value) + (reciprocal of ij energy) * (derivative)
            dv2dx(i) = dv2dx(i) + v2(j)*dv2dxij(nn)
            dv2dy(i) = dv2dy(i) + v2(j)*dv2dyij(nn)
            dv2dz(i) = dv2dz(i) + v2(j)*dv2dzij(nn)

         end do
      end do
C So now we have all derivative totals.
C Lets sum up and return these values for further calculation

C (Original Comment)      
c     calculate total derivatives

      do i=1,natoms
C     Overall atoms
         g(i) = 2.d0 * dv1dx(i) + dv2dx(i) 		! Gradient of x
         g(i+natoms) = 2.d0 * dv1dy(i) + dv2dy(i)	! Gradient of y
         g(i+2*natoms) = 2.d0 * dv1dz(i) + dv2dz(i)	! Gradient of z
      end do

C     Set total energy
      f = v

      return

C     Error message as cluster exploded
 99   iflag = .true.

      return 
      end
