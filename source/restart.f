c    Genetic Algorithm Program
c    Develoepd by the Birmingham University Cluster Group, Prof.
c    Roy Johnston
c    Commenting added by Andrew Logsdail, Oct 2010

C    This is our subroutine to create the DFT restart files. SH 01/11
C    UPDATES 7/11 AJL: - Now just saves current population,
C                        not offspring and parents too.
C    UPDATES 7/11 AJL: - Moved read restart into this file
C                      - renamed restart_DFT to write_restart

          SUBROUTINE write_restart(pop,q)
          use commons, only: natoms,nclust,nmax,ncmax
          implicit none
C          integer noff             ! Number of offspring
C          integer nmut             ! Number of mutants
          integer k
          integer l

          real*8 pop((3*nmax+1)*ncmax)       ! Offspring
          real*8 q(nmax*ncmax)               ! Charge (element identifier)

          logical fexist
 
          inquire(file='popRestart',exist=fexist)

          if (fexist) then
             open(unit=28,file='popRestart',status='old')
             close(28,status='delete')
          end if   

          open(unit=28,file='popRestart',status='new')
C          do l = 1,(nclust+noff+nmut)
           do l = 1, nclust

            write(28,'(f18.10)')pop(l*(3*natoms+1))
            do k = 1,natoms
              write(28,'(f14.10,2x,f14.10,2x,f14.10,2x,f14.10)')
     &              q(k+(l-1)*natoms),
     &              pop(k+(l-1)*(3*natoms+1)),
     &              pop(k+natoms+(l-1)*(3*natoms+1)),
     &              pop(k+2*natoms+(l-1)*(3*natoms+1))
            end do
          end do
          close(28,status='keep')      

          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         SUBROUTINE read_restart(nclust,natoms,nstride,pop,q,energy)

      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
         implicit none
C         include "ga.inc"

         integer natoms           ! Number of atoms
         integer nclust           ! Number of clusters
         integer k
         integer l
         integer nstride	  ! Stride length between clusters

         real*8 pop((3*nmax+1)*ncmax)      ! Offspring
         real*8 q(nmax*ncmax)              ! Charge (element identifier)
         real*8 energy(ncmax)		   ! Energy of clusters

         logical fexist

         inquire(file='popRestart',exist=fexist)

         if (fexist) then
            write(*,*) 'Reading coordinates from restart file'

            open(unit=77,file='popRestart',status='old')

C At this stage we only want number of clusters, so we'll just read/write that infomation
            do l = 1,nclust

               read(77,'(f18.10)')pop(l*nstride)
               do k = 1,natoms
                  read(77,'(f14.10,2x,f14.10,2x,f14.10,2x,f14.10)')
     &            q(k+(l-1)*natoms),
     &            pop(k+(l-1)*nstride),
     &            pop(k+natoms+(l-1)*nstride),
     &            pop(k+2*natoms+(l-1)*nstride)
               end do
            end do
         
            close(77,status='keep')

            do l = 1,nclust
               energy(l) = pop(l*nstride)
C                  PRINT *, energy(l)
            end do

         else 
            write(*,*) 'Error reading from popRestart'
            write(*,*) 'File does not exist'
         end if

         END
