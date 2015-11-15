c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     (Original Comments)
c     contains fitness calculating routines :
c     exp_fitness
c     cos_fitness
c     lin_fitness
c     tanh_fitness
c     pow_fitness
c     pred_fit - gives 'global minimum' value zero fitness

C Subroutine exp_fitness

      SUBROUTINE exp_fitness(nclust,fit,energy)
      implicit none

C Inherited Integer Values
      integer nclust
C Internal Integer Values
      integer i
C Inherited Real Values
      real*8 fit(nclust)
      real*8 energy(nclust)
C Internal Real Values
      real*8 best
      real*8 erange
         
C (Original Comment)
c     Calculates fitness
C Energies are sorted at this point

      best = energy(1)			! Best energy
      erange = energy(nclust)-energy(1)	! Energy range

      do i=1,nclust
C For all clusters exponential fitness
C
C fitness = e^(3*(best-current/range))

         fit(i) = dexp(3.d0*((best-energy(i))/erange)) 
          
      end do

      return	! Return
      end

C Subroutine cos_fitness

      subroutine cos_fitness(nclust,fit,energy)
      implicit none

C Inherited Integer Values
      integer nclust
C Internal Integer Values
      integer i
C Inherited Real Values
      real*8 fit(nclust)
      real*8 energy(nclust)
C Internal Real Values
      real*8 best
      real*8 erange
      real*8 cosx        

C (Original Comment) 
c     Calculates fitness
C Energies are sorted at this point

      best = energy(1)                  ! Best energy
      erange = energy(nclust)-energy(1) ! Energy range

      do i=1,nclust
C For all clusters cos fitness
C
c fitness = (cos[((best-current)/range)*1.4])^2

         cosx = cos(((best-energy(i))/erange)*1.4)
         fit(i) = cosx*cosx 
C Debug
c         print *,i,fit(i),energy(i)
          
      end do

      return		! Return
      end
      
C Subroutine lin_fitness

      SUBROUTINE lin_fitness(nclust,fit,energy)
      implicit none

C Inherited Integer Values
      integer nclust
C Internal Integer Values
      integer i
C Inherited Real Values
      real*8 fit(nclust)
      real*8 energy(nclust)
C Internal Real Values
      real*8 best
      real*8 erange

C (Original Comment) 
c     Calculates fitness
C Energies are sorted at this point

      best = energy(1)                  ! Best energy
      erange = energy(nclust)-energy(1) ! Energy range

      do i=1,nclust
C For all clusters linear fitness
C
c fitness = 1 - 0.7*((current-best)/range)

         fit(i) = 1 - 0.7d0*((energy(i)-best)/erange)

      end do

      return	! Return
      end

C     Subroutine tan_fitness

      SUBROUTINE tanh_fitness(nclust,fit,energy)
      implicit none

C Inherited Integer Values
      integer nclust
C Internal Integer Values
      integer i 
C Inherited Real Values
      real*8 fit(nclust)
      real*8 energy(nclust)
C Internal Real Values
      real*8 best
      real*8 erange
      
C (Original Comment) 
c     Calculates fitness
C Energies are sorted at this point

      best = energy(1)                  ! Best energy
      erange = energy(nclust)-energy(1) ! Energy range

      do i=1,nclust
C For all clusters tan fitness
C
c fitness = 0.5 * (1 - tan(2*((current-best)/range)-1

         fit(i)=0.5d0*(1-dtanh(2.d0*((energy(i)-best)/erange)-1.d0))

      end do

      return	! Return
      end

C Subroutine pow_fitness

      SUBROUTINE pow_fitness(nclust,fit,energy)
      implicit none
      
C Inherited Integer Values
      integer nclust
C Internal Integer Values
      integer i
C Inherited Real Values
      real*8 fit(nclust)
      real*8 energy(nclust)
C Internal Real Values
      real*8 best
      real*8 erange

C (Original Comment) 
c     Calculates fitness
C Energies are sorted at this point
      
      best = energy(1)                  ! Best energy
      erange = energy(nclust)-energy(1) ! Energy range
      
      do i=1,nclust
C For all clusters power fitness
C
c fitness = 1 - ((current-best)/range)^2

         fit(i)= 1 - ((energy(i)-best)/erange)**2.d0

      end do

      return	! Return
      end

C Subroutine pred_fit
C gives 'global minimum' value zero fitness

      subroutine pred_fit(nclust,fit,energy,gmin)
      implicit none
     
C Inherited Integer values 
      integer nclust
C Internal Integer values
      integer i
      integer ic
C Inherited real values
      real*8 fit(nclust)
      real*8 energy(nclust)
      real*8 gmin
C Internal real value      
      real*8 best
      real*8 erange

C Set starting values
      ic = 0
      i = 1

      do while (ic .eq. 0)
C Whilst ic = 0
         if (energy(i) .gt. (gmin+1.e-6)) then
C If energy is greater than max
            best = energy(i)	! Store current energy
            ic = 1		! Note high energy
            print *,i,best,energy(i)	! Print to screen
         else 
            i = i+1	! Keep going til we reach the end
         end if
      end do

C Energy range is last - best
      erange = energy(nclust) - best

      do i=1,nclust
C For all clusters
         if (energy(i) .gt. (gmin+1.e-6)) then
C If energy is greater than maximum
C Set fitness of this current value as exponential value
            fit(i) = dexp(3.d0*((best-energy(i))/erange)) 
         else
C Otherwise fitness is 0 i.e. won't be selected again
            fit(i) = 0.d0
         end if
      end do

      do i=1,nclust
C For all clusters, print energy and minimum energy
         print *,energy(i),fit(i),gmin
      end do

      return
      end

