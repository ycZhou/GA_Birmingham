c This is a program for calculating moments of inertia for n atom clusters. It provides following steps:
c -Transformation of coordinates to centre of mass coordinates
c -Generation of the inertia matrix in centre of mass coordinates
c -Diagonalization of the inertia matrix to obtain moments of inertia along the principal axis.   

       SUBROUTINE moi(mass_a, mass_b, natoms, nstride, na, nb,
     &        pop, q, nclust, inertia)

       IMPLICIT NONE

c      .. Input ..
       DOUBLE PRECISION  pop(*)        !geometry matrix input
       DOUBLE PRECISION  q(*)          !charge matrix input
       DOUBLE PRECISION  mass_a        !mass of atom type a in amu
       DOUBLE PRECISION  mass_b        !mass of atom type b in amu

       INTEGER           natoms        !number of atoms 
       INTEGER           nclust        !number of cluster
       INTEGER           nstride       !3*natoms+1
       INTEGER           na            !number of atoms type a
       INTEGER           nb            !number of atoms type b

c      .. Arguments ..
       DOUBLE PRECISION   A(3,3)       !moment of inertia matrix
       DOUBLE PRECISION   P(nstride)   !input from POP, x1..n,y1..n,z1..n,E,x2..n,y2..n,... 
       DOUBLE PRECISION   sum
       REAL*8             I_xx, I_xy, I_xz, I_yy, I_yz, I_zz !matrix elements of A(3,3)
       REAL*8             sum_mx, sum_my, sum_mz !sums for moment of inertia calculations
       DOUBLE PRECISION   QX(natoms)   !charge array, discrimination of different atom types
       DOUBLE PRECISION   S(3,3)       !matrix employed during matrix diagonalization 
       DOUBLE PRECISION   R(3)         !center of mass vector
       DOUBLE PRECISION   Mass         !mass of atom type b in amu
       DOUBLE PRECISION   temp         !temporary variable 
       DOUBLE PRECISION   W(3)         !eigenvalue of the moment of inertia matrix A

c     .. Integers ..
      INTEGER             i, j, n      !running indices

c     .. Output ..
      DOUBLE PRECISION   inertia(3*nclust)

!      na = 2
!      nb = 1
!      natoms = 3
!      nstride = 3*natoms+1


      do 10 n = 1, nclust

      QX = 0.d0
       P = 0.d0

c writing the QX matrix from input q

      do 20 i = (n-1)*natoms+1, n*natoms 
       j = i - (n-1)*natoms 
       QX(j) = q(i)
 20   end do

!      QX(1) = 2.0
!      QX(2) = 1.0
!      QX(3) = 1.0
 
c writing the P matrix from input pop

      j = 0      
  
      do 30 i = (n-1)*nstride+1, n*nstride
       j = i - (n-1)*nstride
       P(j) = pop(i)
 30   end do

!      P(1) = 0.0
!      P(2) = 0.0
!      P(3) = 0.0
!      P(4) = 0.0
!      P(5) = 0.0
!      P(6) = 1.27

!      P(1) = 0.0
!      P(2) = 0.0
!      P(3) = 0.0
!      P(4) = 0.0
!      P(5) = -0.7575
!      P(6) = 0.7575
!      P(7) = 0.0652
!      P(8) = -0.5213
!      P(9) = -0.5212

c set values for Mass and R

          Mass = na * mass_a + nb * mass_b
          R = 0.d0

c calculate the center of mass vector R

      do i = 1, natoms
       if (QX(i) .eq. 1.d0) then
        R(1) = R(1)+P(i)*mass_a/Mass
        R(2) = R(2)+P(i+natoms)*mass_a/Mass
        R(3) = R(3)+P(i+2*natoms)*mass_a/Mass
       else
        R(1) = R(1)+P(i)*mass_b/Mass
        R(2) = R(2)+P(i+natoms)*mass_b/Mass
        R(3) = R(3)+P(i+2*natoms)*mass_b/Mass
       end if
      end do

c Convert to center of mass coordinate system

      do i = 1, natoms
       P(i) = (P(i)-R(1))
       P(i+natoms) = (P(i+natoms)-R(2))
       P(i+2*natoms) = (P(i+2*natoms)-R(3))
      end do

c     sum over m*x²
 
      sum = 0.d0
       do 40 i = 1, natoms
        if (QX(i) .EQ. 1.0) then
         sum = sum + mass_a*P(i)**2
        else
         sum = sum + mass_b*P(i)**2
        end if
   40  end do
      sum_mx = sum
     
c    sum over m*y²

      sum = 0.d0
       do 50 i = 1, natoms
        if (QX(i) .EQ. 1.0) then
         sum = sum + mass_a*P(i+natoms)**2
        else   
         sum = sum + mass_b*P(i+natoms)**2
        end if
   50  end do
      sum_my = sum

c    sum over m*z²     
      
      sum = 0.d0
       do 60 i = 1, natoms
        if (QX(i) .EQ. 1.0) then
         sum = sum + mass_a*P(i+2*natoms)**2
        else
         sum = sum + mass_b*P(i+2*natoms)**2
        end if
   60  end do
      sum_mz = sum

c     I_xy
   
       sum = 0.d0
       do 70 i = 1, natoms
        if (QX(i) .EQ. 1.0) then
         sum = sum + mass_a*P(i)*P(i+natoms)
        else
         sum = sum + mass_b*P(i)*P(i+natoms)
        end if
   70  end do
      I_xy = -sum

c     I_xz

       sum = 0.d0
       do 80 i = 1, natoms
        if (QX(i) .EQ. 1.0) then
         sum = sum + mass_a*P(i)*P(i+2*natoms)
        else
         sum = sum + mass_b*P(i)*P(i+2*natoms)
        end if
   80  end do
      I_xz = -sum
      
c     I_yz

       sum = 0.d0
       do 90 i = 1, natoms
        if (QX(i) .EQ. 1.0) then
         sum = sum + mass_a*P(i+natoms)*P(i+2*natoms)
        else
         sum = sum + mass_b*P(i+natoms)*P(i+2*natoms)
        end if
   90  end do
      I_yz = -sum

c calculation matrix elements inertia matrix

      I_xx = sum_my + sum_mz
      I_yy = sum_mx + sum_mz
      I_zz = sum_mx + sum_my

c definition of the inertia tensor

      A(1,1) = I_xx
      A(1,2) = I_xy
      A(1,3) = I_xz
      A(2,1) = I_xy
      A(2,2) = I_yy
      A(2,3) = I_yz
      A(3,1) = I_xz
      A(3,2) = I_yz
      A(3,3) = I_zz
 
!      write(*,*) 'A_before', A
      
c diagonalization of 3x3 inertia matrix

      call DSYEVJ3(A, S, W) 

!      write(*,*) 'A', A
!      write(*,*) 'S', S
!      write(*,*) 'W', W

c sort the moments of inertia 

      do 100 i = 1,3
       do 110 j = i,3
        if (W(i) .gt. W(j)) then
          temp = W(i)
          W(i) = W(j)
          W(j) = temp
        end if
 110   end do
 100  end do

c calculate relative moments of inertia

      inertia(2*nclust+n) = W(3)/W(1)
      inertia(nclust+n) = W(2)/W(1)
      inertia(n) = W(1)/W(1)

!      write(*,*) 'moments of inertia', W 

   10  end do 

c      write(*,*) 'inertia', inertia

c end of moi

      end 
