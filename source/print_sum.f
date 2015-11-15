c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     (Original Comments)
c     routines for ouput :
c     print_sum - summary of a generation
c     print_init - writes out xyz file of every cluster in initial pop
c     print_current - writes out lowest energy cluster in pop, if lower than 
c     previous
c     av - calculates average energy of population

C Subroutine print_sum
C summary of a generation

      subroutine print_sum(nclust,energy,i)
      implicit none

C Inherited integer values
      integer nclust	! Number of clusters
      integer i		! Generation Number
C Inherited real values
      real*8 energy(nclust)	! Energies
C Internal real values
      real*8 av
      real*8 avg
C Import averaging function
      external av
      
C Print to screen values
      write(*,'('' Generation : '',i3)')i	  		! Generation Number
      write(*,'('' Lowest Energy  = '',f13.6)')energy(1)   	! Best
      write(*,'('' Highest Energy = '',f13.6)')energy(nclust)	! Worst
      avg = av(nclust,energy)
      write(*,'('' Average Energy = '',f13.6,/)')avg		! Average
      
      write(4,'(1x,i3,3x,f13.6,3x,f13.6,3x,f13.6)')
     & 	   i,energy(1),energy(nclust),avg			! Write to file

      return	! Return
      end

C Subroutine print_sum_ne
C Same as above, but with supplied minimum

      subroutine print_sum_ne(nclust,energy,i,best)
      implicit none
            
C Inherited integer values
      integer nclust	! Number of clusters
      integer i		! Generation number
C Inherited real values
      real*8 energy(nclust)
      real*8 best
C Internal real values
      real*8 av
      real*8 avg
C Import averaging function
      external av

C Print to screen values
      write(*,'('' Generation : '',i3)')i	  		! Generation Number
      write(*,'('' Lowest Energy  = '',f13.6)')best		! Best
      write(*,'('' Highest Energy = '',f13.6)')energy(nclust)	! Worst
      avg = av(nclust,energy)
      write(*,'('' Average Energy = '',f13.6,/)')avg		! Average
      
      write(4,'(1x,i3,3x,f13.6,3x,f13.6,3x,f13.6)')
     & 	   i,best,energy(nclust),avg				! Write to file

      return	 ! Return
      end

C Averaging function, used by above      

      function av(nclust,energy)
      implicit none

C Inherited integer values
      integer nclust	! Number of clusters
C Internal integer value
      integer i
C Inherited real values
      real*8 energy(nclust)	! Energies
C Internal real value
      real*8 total
      real*8 av

C (Original Comments)
c     calculates the average energy of the population
      
      total = 0.0d0

      do i=1,nclust
C Sum totals of clusters
         total = total + energy(i)
      end do

      av = total/dble(nclust)		! Divide total by number of clusters

      return		! Return average value
      end

C Subroutine write_gen
C Writes generation to file? 
      
      subroutine write_gen(natoms,nclust,nstride,pfunc,nels,pop,energy,
     &           q,el,namea,nameb)
      implicit none

C Inherited integer values
      integer natoms	! Number of atoms
      integer nclust	! Number of clusters
      integer nstride	! Step size
      integer pfunc	! Potential
      integer nels	! Number of elements
C Internal integer values
      integer i
      integer j
C Inherited real values
      real*8 pop(nstride*nclust)	! Coordinates
      real*8 q(natoms*nclust)		! Charges
      real*8 energy(nclust)		! Energies
C Inherited string values
      character*(*)el
      character*(*)namea
      character*(*)nameb

      if ( nels .eq. 2 ) then
C (Original Comments)
c     two element clusters

         if ( pfunc .eq. 3 .or. pfunc .eq. 5) then
C (Original Comments)
c     ionic clusters

            do i=1,nclust
C For all clusters

               write(7,'(f13.6)')energy(i)	! Print energy

               do j=1,natoms
C For all atoms in this cluster
               
                  if (q(j+(i-1)*natoms) .gt. 0.d0) then
C If charge is greater than 0 type A
                     write(7,'(a,1x,3f13.8)')namea,pop(j+(i-1)*nstride),
     &                    pop(j+natoms+(i-1)*nstride),
     &                    pop(j+2*natoms+(i-1)*nstride)
                  else
C Charge is less than 0 sowe have type B
                     write(7,'(a,1x,3f13.8)')nameb,pop(j+(i-1)*nstride),
     &                    pop(j+natoms+(i-1)*nstride),
     &                    pop(j+2*natoms+(i-1)*nstride)
                  end if
                  
               end do
            end do

         else if ( pfunc .eq. 4 .or. pfunc .eq. 6 ) then
C We have Gupta Potential Type with bimetallics
            do i=1,nclust
C For all clusters

               write(7,'(f13.6)')energy(i)	! Print energy

C (Original comments)
c     alloy clusters (gupta potential )               

               do j=1,natoms
C For all atoms in cluster
                     
                  if (q(j+(i-1)*natoms) .eq. 1.d0) then
C If charge is 1 we have type A
                     write(7,'(a,1x,3f13.8)')namea,pop(j+(i-1)*nstride),
     &                    pop(j+natoms+(i-1)*nstride),
     &                    pop(j+2*natoms+(i-1)*nstride)
                  else
C Else we have type B
                     write(7,'(a,1x,3f13.8)')nameb,pop(j+(i-1)*nstride),
     &                    pop(j+natoms+(i-1)*nstride),
     &                    pop(j+2*natoms+(i-1)*nstride)
                  end if
                  
               end do
            end do

         end if

      else
C (Original Comments)
c     single element clusters

         if ( pfunc .eq. 4 .or. pfunc .eq. 6 ) then
C Gupta potential for monometallics
C (Original Comment)
c     slight dodgyness with element arrays (namea not el)

            do i=1,nclust
C For all clusters
               write(7,'(f13.6)')energy(i)
C Print energies               
               do j=1,natoms
C Print atom coordinates                     
                  write(7,'(a,1x,3f13.8)')namea,pop(j+(i-1)*nstride),
     &                 pop(j+natoms+(i-1)*nstride),
     &                 pop(j+2*natoms+(i-1)*nstride)

               end do
            end do

         else
C For all other cases
            do i=1,nclust
C Print cluster energy values
               write(7,'(f13.6)')energy(i)
               
               do j=1,natoms
C Print atom coordinates                     
                  write(7,'(a,1x,3f13.8)')el,pop(j+(i-1)*nstride),
     &                 pop(j+natoms+(i-1)*nstride),
     &                 pop(j+2*natoms+(i-1)*nstride)
                  
               end do
            end do

         end if
      end if

      return	! Return

      end
            
C Subroutine write_en
C Print energies to screen
           
      subroutine write_en(nclust,n,energy) 
      implicit none

C Inherited integer values
      integer nclust	! Array of energies
      integer n		! Current generation
C Internal integer values
      integer i
C Inherited real values
      real*8 energy(nclust)	! Energies
      
      write(8,'(i3,100f13.6)') n,(energy(i), i=1,nclust)	! Print all energies to screen

      return	! Return
      end
    
