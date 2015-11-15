c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C This is our minimising subroutine, calling the LBFGS

      SUBROUTINE minimiser(natoms,na,nb,idum,pfunc,idamp,ihard,nels,off,
     &           qoff,alpha,r_ab,qa,qb,B_ab,rho_ab,B_bb,rho_bb,de,re,a2,
     &           a3,coeff,rcut2,range,nspec,aij,vij,pij,qij,r0,
     &           c1_sisi,c2_sisi,c3_sisi,c4_sisi,
     &           c1_sio,c2_sio,c3_sio,c4_sio,
     &           c1_oo,c2_oo,c3_oo,c4_oo,
     &           ru_sisi,ru_sio,ru_oo,cutoff_gupta,cutoff_start,
     &           cutoff_end,a3a,a4a,a5a,x3x,x4x,x5x,j,k,
     &           namea,nameb,mass_a,mass_b)

      use hectorstuff
      use DFT, only: dft_type, submit
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
      implicit none
C Setup parameter (Limited variable)
      integer mmax
      parameter(mmax=20)

C Integer Values
      integer natoms		! Number of atoms	
      integer na		! Number of atoms a
      integer nb		! Number of atoms b
      integer idum		! Seed for random number generator
      integer pfunc		! Function Type
      integer idamp
      integer ihard
      integer nspec		! Number of species
      integer nels		! Number of elements
C      integer grid      ! grids for GPAW; to be multiples of 8
C Minimiser Parameters
      integer m
      integer iprint
      integer i
      integer j                 ! is 'i' in ga.f and calls DFT
      integer k                 ! is 'j' in ga.f and calls DFT
      integer nbd(3*nmax)	! Bounds for minimiser
      integer isave(44)
      integer iwa(9*nmax)

C Real Values
      real*8 off(3*natoms+1)	! Offpsring
      real*8 qoff(natoms)	! Charge on offspring
C Morse Potential Parameter
      real*8 alpha
C Unknown parameters
      real*8 r_ab
      real*8 qa			! Charge on atom A
      real*8 qb			! Charge on atom B
      real*8 B_ab
      real*8 rho_ab
      real*8 B_bb
      real*8 rho_bb
C MM Parameters
      real*8 de
      real*8 re
      real*8 a2
      real*8 a3
      real*8 coeff(11)
      real*8 rcut2
      real*8 range	! Atom Range
C Gupta Parameters (Arrays)
C These are arrays as they are configured for bimetallics
C Thus holding interactions A-A, A-B, B-A and B-B
      real*8 aij(nspecmax2)
      real*8 vij(nspecmax2)
      real*8 pij(nspecmax2)
      real*8 qij(nspecmax2)
      real*8 r0(nspecmax2)
C Si Parameters
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
      real*8 ru_sio
      real*8 ru_oo 
C Cutoff Parameters
      real*8 cutoff_start
      real*8 cutoff_end
      real*8 a3a	! Repulsive polynomial coefficients
      real*8 a4a
      real*8 a5a
      real*8 x3x	! Attractive polynomial coefficients
      real*8 x4x
      real*8 x5x
C Minimiser Parameters
      real*8 f			! Energy value return from calculator
      real*8 factr
      real*8 pgtol		! Minimum residual gradient for convergence
      real*8 g(3*nmax)		! Gradients
      real*8 l(3*nmax)
      real*8 u(3*nmax)
      real*8 dsave(29)
      real*8 wa(2*mmax*3*nmax+12*nmax+12*mmax*mmax+12*mmax)
      real*8 ran3		! Random Number Generator
      real*8 mass_a ! Mass of element a -> write_DFT.f
      real*8 mass_b ! Mass of element b -> write_DFT.f
      real*8 size_cell ! Cell size for DFT call
      real*8 jobDone   ! Is DFT job done ?
C      real*8 uc         ! box size for GPAW
c      integer nproc     !number of procs (HECTOR ONLY)
c      integer ntask     !number of tasks (HECTOR ONLY)


C Boolean operator
C For cutoff
      logical cutoff_gupta
C For minimiser
      logical lsave(4)
      logical iflag
      logical ok

C External parameter
      external ran3	! Random numbers

C String parameters
C For minimiser
      character*60 task		 ! Feedback from minimiser
      character*60 csave         
      character*60 prefix        ! Prefix of filename for DFT call
      character*(*) namea        ! Name of element a
      character*(*) nameb        ! Name of element b
      character*4  NPROCSTR      !string variable for number of procs (HECTOR)
      character*4  NTASKSTR      !string variable for number of ppn (HECTOR)




CCCCCCCCCCCCCCCCCCCCC Variables all declared CCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCC Now to initialise      CCCCCCCCCCCCCCCCC

C     (Original Comments)
c     set up minimisation options

      m=5
      factr=1.0d+7
      pgtol=1.0d-5		! Residual gradient
      iprint = -1

C     (Original Comment)
c     no bounds
 
      do i=1,3*natoms
         nbd(i) = 0
      end do
   
C     Setup to run LBFGS-B Routine
      task = 'START'

C     (Original Comments)
c     start  L-BFGS-B Routine

      if (pfunc .eq. 1 ) then

C     (Original Comments)
c     morse potential
C     Minimiser whilst not minimised
 
         do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)

C           Calculate forces         
            if (task(1:2) .eq. 'FG') call morse(natoms,alpha,
     &           off,f,g)
         
	 end do
               
      else if ( pfunc .eq. 2 ) then

C     (Original Comments)
c     murrell-mottram
C     Minimiser while not minimised

         do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)
            
C           Calculate forces
            if (task(1:2) .eq. 'FG') call murrmott(natoms,idamp,
     &           ihard,de,re,a2,a3,coeff,rcut2,off,f,g)

         end do
         
      else if ( pfunc .eq. 3 ) then

C     (Original Comments)
c     rigid-ion
C     Minimiser while not minimised
         
 99      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)
         
C           Calculate forces
            if (task(1:2) .eq. 'FG') call ionic(natoms,B_ab,rho_ab,
     &           B_bb,rho_bb,off,f,g,qoff,qa,qb)

         end do

C     (Original Comment)
c     generate new cluster if collapsed - dodgy potential
C Checks if forces are less than arbitrary value.
C If so recreate atoms using ionic method (gen_pop_ionic)

         if (f .lt. -40.d0*dble(natoms)/2.d0) then
            write(*,
     &         '('' Cluster has collapsed, generating replacement'',/)')

C For every coordinate in the system
            do i=1,3*natoms
C Generate new offspring coordinates
               off(i) = r_ab*ran3(idum)*range
            end do
C For all atoms A
            do i=1,na
C Set charges of A
               qoff(i) = qa
            end do
C And for all of atoms B
            do i =1,nb
C Set charges of B
               qoff(i+na) = qb
            end do

C Set task to start again, and return to rigid-ion minimiser
            task = 'START'
            goto 99
         end if
         
      else if ( pfunc .eq. 4 ) then

C     (Original comment)
c     gupta potential
C     Minimiser until minimised
            
 111        do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')

            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)

C           Calculate forces
            if (task(1:2) .eq. 'FG') call gupta(natoms,nspec,off,qoff,f,
     &           g,aij,vij,pij,qij,r0,iflag,cutoff_gupta,cutoff_start,
     &           cutoff_end,a3a,a4a,a5a,x3x,x4x,x5x)

C Check if cluster has exploded
C If so generate a new cluster
            if (iflag) then
               write(*,
     &       '('' Cluster has exploded, generating replacement'',/)')
               
C If we are dealing with a bimetallic cluster...
               if ( nels .ge. 2 ) then
C For all atoms assign random coordinates
                  do i=1,3*natoms
                     off(i) = ran3(idum)*range
                  end do
C Then for atoms a assign charges of a
                  do i=1,na
                     qoff(i) = qa
                  end do
C and for atoms b assign charge of b
                  do i =1,nb
                     qoff(i+na) = qb
                  end do
C And restart
                  task = 'START'
                  goto 111
               else
C Else assign random coordinates for all
                  do i=1,3*natoms
                     off(i) = ran3(idum)*range
                  end do
C And restart
                  task = 'START'
                  goto 111
               end if

            end if

         end do
            
C Print minimisation to save file
c     print *,isave(34)

      else if ( pfunc .eq. 5 ) then

C     (Original Comment)
c     TTAM potential for Si-O
C     Minimiser until minimised

199      do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)
         
C           Calculate forces
            if (task(1:2) .eq. 'FG') 
     &       call silica(natoms,na,nb,c1_sisi,c2_sisi,c3_sisi,c4_sisi,
     &       c1_sio,c2_sio,c3_sio,c4_sio,c1_oo,c2_oo,c3_oo,c4_oo,
     &       ru_sisi,ru_oo,ru_sio,off,f,g,qoff,qa,qb,ok)

C     (Original comment)
c     generate new cluster if unphysical                   
           if (.not.ok) then 
C             write(*,'('' Unphysical cluster, generating replacement'',/)')
              do 
C For all atoms generate new coordinates in box
                do i=1,3*natoms
                   off(i) = 3*ran3(idum)*range
                end do
C Assign charges to atoms a
                do i=1,na
                   qoff(i) = qa
                end do
C And assign charges to atoms b
                do i =1,nb
                   qoff(i+na) = qb
                end do
                ok=.true.

C Check if structures are acceptable
                call check_unphys_silica(ok,na,nb,off,
     &                qoff,qa,qb,ru_sisi,ru_sio,ru_oo)
              if (ok) exit 
              end do

C And restart
             task = 'START'
             goto 199
           end if
         end do

      else if ( pfunc .eq. 6 ) then

! SPACE FOR DFT CALLS
         jobDone = 0.d0
         do while (1.0d0 .gt. jobdone)

C Create input file    
            call write_DFT(natoms,na,nb,off,qoff,j,k,
     &               namea,nameb,
     &               mass_a,mass_b,nspec,
     &               prefix,size_cell,
     &               range,r_ab)
       
C System calls to run DFT codes
            if (dft_type .eq. 1) then
            !! QE RUN COMMAND
            CALL system(trim(submit)//' < '//trim(prefix)//'.in'
     &           //' > '// trim(prefix)//'.out')


            else if (dft_type .eq. 2) then

            !! GPAW RUN COMMAND
            call system('ulimit -t 10800;'//
     &      ' mpirun -np `wc -l < $PBS_NODEFILE`'//
     &      ' gpaw-python '//trim(prefix)//'.in > gpaw.log')
C     &      ' python '//trim(prefix)//'.in > gpaw.log') ! Used for testing


            else if (dft_type .eq. 3) then

            !! G09 RUN COMMAND
C         cmd = 'ulimit -t 18000; g09 <' // trim(name1) // '.in > '
C     &         // trim(name1) // '.log'
CC        print *, cmd
C         call system(cmd)

            call system('ulimit -t 18000;'//
     &     ' g09 <'//trim(prefix)//'.in > '//trim(prefix)//'.out')


            else if (dft_type .eq. 4) then
            !! NW RUN COMMAND
             call system(trim(submit)//' '//trim(prefix)//'.nw'//
     &       ' > '//trim(prefix)//'.out')


            !! VASP RUN COMMAND
          else if (dft_type .eq. 5) then
            call system('cd '//trim(prefix)//'/ && '//trim(submit))
          end if

C Gather information from output file
            call output_DFT(natoms,prefix,qoff,off,namea,nameb,
     &            size_cell,iflag,jobDone)

! If not converged: (1) Generate random strucutre.
! Resubmit new structure and start again!
            if (iflag) then
               write(*,'('' Calculation not converged!'',/)')
               jobDone = 0.d0
! Choose Path
               write(*,'('' Resubmit new Geometry!'',/)')

               do i=1,3*natoms
                  off(i) = ran3(idum)*range*r_ab
               end do

               if ( nels .ge. 2 ) then
                  do i=1,na
                     qoff(i) = qa
                  end do
                  if (nb .gt. 0) then
                     do i=1,nb
                        qoff(i+na) = qb
                     end do
                  end if
               end if
            end if
! End of resubmit new structure appenditure
         end do
      end if

C Assign energy to the end offspring?
C Unsure on this one!

      if (pfunc .ne. 6) off(3*natoms+1) = f
      
      return
      end 
