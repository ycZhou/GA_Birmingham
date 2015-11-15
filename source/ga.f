c     Genetic Algorithm Program
c     Developed by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C For full reference please see Dalton Trans article (2003, R. L. Johnston)

C     Note to self: Need to figure OMP implementation - at present causes a segmentation faults at OMP_SET_NUM_THREADS()

      PROGRAM Genetic_Algorithm

c     Main program executable

      use commons
      use potentials
      implicit none

c     Integer values for GA
      integer nmut	! Number of Mutations
      integer ifilelen

C     Integer values for internal calculation
      integer i		!
      integer j		!
      integer k		!
      integer l		!
      integer iterm	!
      integer mflag	!
      integer plen	!
      integer str(2)	!
      integer narg	! Number of command line arguments
      integer ilen	! Length of String variable
      integer iargc	!
      integer mat	!
      integer mut	!
      integer nn	!
      integer nmin	!
C      integer grid      ! grids for GPAW; to be multiples of 8

C     Real Numbers
      real*8 mfact	!
      real*8 rprev	!
      real*8 rcurr	!
      real*8 eold	!
      real*8 ran3	! Random number
      real*8 oldCluster ! Lowest energy in old generation
C      real*8 uc         ! box size for GPAW (dft_type = 2)
      real*8 inertia(1000)
C     Boolean values
      logical gflag	!
      logical ok	!
      logical fexist
      
C     Strings and / or Character Representations
      character*20 potential 	! Potential file name
      character*1 rep	!
      character*10 pot	!
      character*20 ifile	! Input file name
      character*20 prog		!

C ADDED BY ANDREW LOGSDAIL. THIS FILE NEEDS SOME SERIOUS TIDYING!!
C Additions added December 2008 to incorporate cutoff in calculations
! MTO: All of these arrays were originally length nspecmax2.
! MTO: For now they're all fixed-length.
      real*8 cutoff_start(9)	! Array of real numbers for cutoff starting distance
      real*8 cutoff_end(9)	! Array of real numbers for cutoff ending distance
      real*8 a3a(9)		! Coefficient a3
      real*8 a4a(9)		! Coefficient a4
      real*8 a5a(9)		! Coefficient a5
      real*8 x3x(9)		! Coefficient x3
      real*8 x4x(9)		! Coefficient x4
      real*8 x5x(9)		! Coefficient x5
      
C     External variables
      external ran3	! Random number generator
      integer lstr	! Length of string
      external OMP_SET_NUM_THREADS	! Set number of OMP threads
     
      real*8 temp(20) 
C     Definition of Pi
      pi = 4.d0*atan(1.d0)
C     Definition of Processors
      processors = 1      

c     (Original comment)
c     read data from input file given on command line

      narg = iargc()

C     Respond if number of arguemnts is incorrect from program

      if (narg .ne. 1) then
         call getarg(0,prog)
         ilen = lstr(prog)						! Get length of input string
         write(*,'('' Usage : '',a,'' <input file>'')')prog(1:ilen)	! Print error
         stop								! Terminate
      end if

C     If correct, get input filename from input arguments

      call getarg(1,ifile)	! Assign input file to ifile
      ilen = lstr(ifile)	! Get length of input file name

C     (Original Comment)
c     open file on unit 2 (expected in routine tokeniser)

C     Open input file retrieved from nargs
      open(unit=2,status='old',file=ifile(1:ilen),iostat=narg)

C     Throw error if no arguments provided
      if (narg .ne. 0) then
         write(*,'(/,''ERROR, cannot open file : '',a,/)')ifile(1:ilen)	! Error message
         stop								! Terminate
      end if

C     Else read input file printed to screen....
      write(*,'(/,''Read input file : '',a,/)')ifile(1:ilen)		! Read file

c     run ends here if errors found in input
C     ** Read in variables from input file **
C     See lib/source/inter.f for further detail
 
      call interpret(potential)
      
      close(unit=2)	! Close input file

C Disabled due to segmentation fault?
C No solutions in my head, or on google!
C Set number of processors
C       write(*,*)processors
C      if (processors .gt. 1) then 
C        call OMP_SET_NUM_THREADS(processors)
C      end if

C     (Original comment)
c     read parameters for potential
      
      if (ptype .eq. 1) then
         inputfile = 'morse.in'
         ifilelen = 8
         pfunc = 1

C        ** read in morse potential from file **
C        See lib/source/read_morse.f for further detail
         call read_morse(inputfile,ifilelen,alpha,bscale,el)

      else if (ptype .eq. 2) then
C
C        Read in softcoded Murrell Mottram
C
         inputfile = potential 
         ifilelen = lstr(potential)
         pfunc = 2

C        ** read in MM potential from file **
C        See lib/source/read_mmparams.f for further detail
         call read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,idamp,
     &        ihard,ram,rcut2,bscale,el,title)

      else if (ptype .eq. 7) then
C
C	 Read in MgO
C
         inputfile = 'MgO.in'
         ifilelen = 6
         pfunc = 3

C        ** read in Rigid-Ion potential from file **
C        See lib/source/read_iparams.f for further detail
         call read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &        B_bb,rho_bb,bscale,namea,nameb)

      else if (ptype .eq. 8) then
C
C        Read in ZnO
C
         inputfile = 'ZnO.in'
         ifilelen = 6
         pfunc = 3

C        ** read in Rigid-Ion potential from file **
C        See lib/source/read_iparams.f for further detail
         call read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &        B_bb,rho_bb,bscale,namea,nameb)

      else if (ptype .eq. 9) then
         inputfile = 'gen_ion.in'
         ifilelen = 6
         pfunc = 3

C        ** read in Rigid-Ion potential from file **
C        See lib/source/read_iparams.f for further detail
         call read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &        B_bb,rho_bb,bscale,namea,nameb)

      else if (ptype .eq. 10) then
C
C        Read in softcoded Gupta
C
         inputfile = potential
         ifilelen = lstr(potential)
         pfunc = 4

C        ** read in Gupta potential from file **
C        See lib/source/read_gupta.f for further detail
         call read_gupta(inputfile,ifilelen,nspec,aij,vij,pij,qij,r0,
     &        bscale,namea,nameb)    

      else if (ptype .eq. 24) then
C
C        Hardcoded SiO
C
         pfunc = 5
         namea   = 'Si'
         nameb   = 'O'
         bscale  = 1.d0
         ru_sisi = 0.0877d0
         ru_sio  = 1.0613d0
         ru_oo   = 1.4760d0

         c1_sisi = 1913.2812d0
         c2_sisi = -537.3124d0
         c3_sisi = 0.06570d0
         c4_sisi = 1.7376d0

         c1_sio  = -956.6406d0
         c2_sio  = -1631.1766d0
         c3_sio  = 0.20851d0
         c4_sio  = 2.9162d0

         c1_oo   = 478.3203d0
         c2_oo   = -4951.9369d0
         c3_oo   = 0.35132d0
         c4_oo   = 4.0948d0

      else if (ptype .eq. 25) then
C
C        Alternative hardcoded SiO
C
         pfunc = 5
         namea   = 'Si'
         nameb   = 'O'
         bscale  = 1.d0
         ru_sisi = 0.977
         ru_sio  = 1.015    
         ru_oo   = 0.985    

         c1_sisi = 1913.2812d0
         c2_sisi = -10303.0160d0
         c3_sisi = 0.201d0
         c4_sisi = 3.2213d0

         c1_sio  = -956.6406d0
         c2_sio  = -1453.9019d0
         c3_sio  = 0.2081d0
         c4_sio  = 2.9043d0

         c1_oo   = 478.3203d0
         c2_oo   = -954.1094d0
         c3_oo   = 0.358d0
         c4_oo   = 4.0918d0

       else if (ptype .eq. 26) then
C
C        Read in DFT parameter file
C
         inputfile = potential 
         ifilelen = lstr(potential)
         pfunc = 6
         bscale = 1.d0 ! Bodge as scaling is not read in for DFT

C        ** read in DFT parameter file **
C        See source/read_DFT.f for further detail
         call read_DFT(inputfile)

      end if

C     range here is the radius of a sphere with volume natoms      
C     this is used for the initial configuration of atoms in on sphere/ in box
      range = dble(natoms)**(1.d0/3.d0)

C     Set realspace gridsizes for DFT calculations if using GPAW (dft_type = 2)
C     Lifted from code donated by Shujiang Yang, ANL
C      if (pfunc .eq. 6 .and. dft_type .eq. 2) then
C This section is added to adjust the GPAW grids to be multiples of 8
C h~0.3; vacuum >= 5.0 on each side (+10.0);
C For CdSe: rab=3.6 (relative); bscale=1.0 not used 
C         uc = range*3.6 + 10.0
C         grid = nint((uc/0.3)/8.0) * 8
C      else
C         uc = 0.d0
C         grid = 0
C      end if
      
c     open summary file
      open(4,status='unknown',file=fname(1:flen)//'.sum')
      
      if ( write_clust ) then
C        Open file to save all clusters
         open(7,status='unknown',file=fname(1:flen)//'.clusters')
         write(7,'(3i3)')natoms,ngen,nclust	! Write out number of atoms, generations and clusters
      end if

      if ( write_energy ) then
C        Open file to save energies
         open(8,status='unknown',file=fname(1:flen)//'.energies')
      end if

      if ( write_stats ) then
C        Open file to save statistics as we go
         open(9,status='unknown',file=fname(1:flen)//'.stats')
      end if

      if (cutoff_gupta) then
C        Read Gupta cutoff parameters
         call read_cutoff(cutoff_start,cutoff_end,a3a,a4a,a5a,
     &        x3x,x4x,x5x)
      end if

C
c     (Original Comments)
c       ** generate initial population **
c     
c     gen_pop_box - puts atoms in box of length natoms**(1/3)
c     gen_pop_sphere - puts atoms on surface of sphere of diameter
c     natoms**(1/3)
C
     
      nstride = 3*natoms+1
      idum = -idum	! This is the seed for the random number generator ran3

C     Check to see if we are restarting the run
C     We'll just code it nicely and put an if statement in...
      
      if (restart) then
C        ** read in old population information from popRestart **
C        See lib/source/restart.f for further detail

            call read_restart(nclust,natoms,nstride,pop,q,energy) 

C            open(unit=77,file='popRestart',status='old')

C At this stage we only want number of clusters, so we'll just read/write that infomation
C             do l = 1,nclust

C               read(77,'(f18.10)')pop(l*nstride)
C               do k = 1,natoms
C                  read(77,'(f14.10,2x,f14.10,2x,f14.10,2x,f14.10)')
C     &            q(k+(l-1)*natoms),
C     &            pop(k+(l-1)*nstride),
C     &            pop(k+natoms+(l-1)*nstride),
C     &            pop(k+2*natoms+(l-1)*nstride)
C               end do
C            end do
C            close(77,status='keep')

C            do l = 1,nclust
C               energy(l) = pop(l*nstride)
C               PRINT *, energy(l)
C            end do

C Set charges
            if (nels .eq. 2) then
               if (pfunc .eq. 4 .or. pfunc .eq. 6) then
c           set qa to 1 and qb to 2, r_ab = 1.d0 for Gupta
                qa = 1.0d0
                qb = 2.0d0
                !! Sven this number needs to be dynamic !!
                !! For DFT it is now read in from the input file !!
                if (pfunc .eq. 4) then
                   r_ab = 1.0d0
                end if
              else if (pfunc .eq. 5) then
                qa = 4.0d0
                qb = -2.0d0
              end if
            end if

      else
C Generate new clusters
         if (nels .eq. 2) then
            if ( pfunc .eq. 3 ) then
C Generate ionic population for Rigid-Ion Potential
               call gen_pop_ionic(na,nb,nclust,nstride,idum,pop,q,qa,qb,
     &              r_ab,range)
c        else if (pfunc .eq. 4) then
            else if (pfunc .eq. 4 .or. pfunc .eq. 6) then
c     set qa to 1 and qb to 2, r_ab = 1.d0 for Gupta
               qa = 1.0d0
               qb = 2.0d0
	       !! Sven this number needs to be dynamic !!
	       !! For DFT it is now read in from the input file !!
               if (pfunc .eq. 4) then 
                  r_ab = 1.0d0
               end if
C Generate Bimetallic population for Gupta Potential (Same method call)
               call gen_pop_ionic(na,nb,nclust,nstride,idum,pop,q,qa,qb,
     &           r_ab,range)
            else if (pfunc .eq. 5) then
               qa = 4.0d0
               qb = -2.0d0
C Generate Population for Silica based Function
              call gen_pop_silica(na,nb,nclust,nstride,idum,pop,q,qa,qb,
     &           ru_sisi,ru_sio,ru_oo,range)
            end if
         else
C Generate Population in a box
            call gen_pop_box(natoms,nclust,nstride,pop,idum,range)
         end if
      
c     relax initial population
      
         nmin = 0
         gflag = .false.

C Parallelisation needs some work. Currently we get segmentation faults
C from the pop() and q() commands, and due to limited stacksize
C
C Can we find a way to allow access to the arrays without issue?
C
C!$OMP PARALLEL DO DEFAULT(PRIVATE) 

C Loop through all clusters
         do i=1,nclust
            j = 0 ! for write DFT call
            write(*,*) 'in ga.f, i (cluster number) is'
            write(*,*)i	! Debug
           write(*,*) 'going into minimiser now'
C And minimise structures
            call minimiser(natoms,na,nb,idum,pfunc,idamp,ihard,nels,
     &      pop((i-1)*nstride+1),q((i-1)*natoms+1),alpha,r_ab,qa,qb,
     &      B_ab,rho_ab,B_bb,rho_bb,de,re,a2,a3,coeff,rcut2,range,
     &      nspec,aij,vij,pij,qij,r0,
     &      c1_sisi,c2_sisi,c3_sisi,c4_sisi,
     &      c1_sio,c2_sio,c3_sio,c4_sio,
     &      c1_oo,c2_oo,c3_oo,c4_oo,
     &      ru_sisi,ru_sio,ru_oo,cutoff_gupta,cutoff_start,
     &      cutoff_end,a3a,a4a,a5a,x3x,x4x,x5x,j,i,
     &      namea,nameb,mass_a,mass_b)
c            write(*,*) 'out of minimiser into ga.f again'
c            write(*,*) 'i is now'
c            write(*,*)i
C Save energy
            energy(i) = pop(i*nstride)

C Number of minimised structures
            nmin = nmin + 1
C Check if we have seen the minimum we are looking for...
            if ( dabs(energy(i)-gmin) .lt. 2.d-6 
     &           .and. .not. gflag) then
                 write(*,'(a,1x,i4)')
     &           'GMIN found by minimization number ',nmin
                 gflag = .true.
            end if
         
         end do 
C!$OMP END PARALLEL DO 
      end if
C (Original Comment)
c     sort energies of initial population
C This subroutine is in pop.f      
     

      call moi(mass_a, mass_b, natoms, nstride, na, nb,
     &  pop, q, nclust, inertia)

 
      call esort(natoms,nstride,pop,q,energy,mutant,nclust,nels)
C Write oldCluster variable for DFT run
      oldCluster = energy(1)
C (Original Comment)      
c     summary of initial population
C This subroutine is in print_sum.f
      
      call print_sum(nclust,energy,0)
      
      if ( write_clust ) then
C (Original Comment)
c     write geometries of current generation to file
C This subroutine is in print_sum.f

         call write_gen(natoms,nclust,nstride,pfunc,nels,pop,energy,
     &        q,el,namea,nameb)

      end if

C Print all energies of population
C This is located in print_sum.f
      if ( write_energy) call write_en(nclust,0,energy)

C (Original Comment)
c     termination counters 
      
      iterm = 0
      rprev = 0.d0
C101   continue
 
C (Original Comment)     
c     begin run of ga
      
      do i=1,ngen
C Write Population information to file for restart for DFT
          if (pfunc .eq. 6) then
           write(*,'('' Writing restart file ...... '',/)')
           call write_restart(pop,q)
          end if

C For number of generations

C         write(*,'('' DEBUG '',/)')
   
         if ( write_stats ) then
C If writing all statistics
            do j=1,nclust+noff
C Set all mutant allocators to zero
               mutant(j) = 0
            end do
         end if
         
         do j=1,nclust
C For all clusters find centre of mass with respect to y coordinate            
            call cofmass(natoms,pop((j-1)*nstride+1))
            
         end do
        
C         write(*,'('' DEBUG '',/)')
 
c     (Original Comment)
c     ** calculates fitness **
C All in fitness.f
         
         if (ftype .eq. 1) then
            call lin_fitness(nclust,fit,energy)		! Linear Fitness
         else if (ftype .eq. 2) then
            call pow_fitness(nclust,fit,energy)		! Power Fitness
         else if (ftype .eq. 3) then
            call exp_fitness(nclust,fit,energy)		! Exponential Fitness
         else if (ftype .eq. 4) then
            call tanh_fitness(nclust,fit,energy)	! Tanh Fitness
         else if (ftype .eq. 5) then
            call cos_fitness(nclust,fit,energy)		! Cos Fitness
         end if
      
c         write(*,'('' DEBUG '',/)')

         do j=1,noff
C For all offspring
c     (Original Comment)
c     ** select a pair of clusters **
C These can all be found in select.f

            if (stype .eq. 1) then
C Roulette Selection
               call roul_select(natoms,nclust,nstride,nels,pop,q,fit,
     &              pair,qp,idum,str)
            else
C Tournament Selection
               call tourn_select(natoms,nclust,nstride,nels,pop,q,
     &              energy,pair,qp,idum,tsize,str)
            end if

C            write(*,'('' DEBUG '',/)')

C (Original Comments)
c rotate clusters - first cluster
            call rotate(natoms,pair,idum,pi)
            
c     note - nstride NOT nstride+1 - cluster is only 3*natoms in pair
C rotate clusters - second cluster
            call rotate(natoms,pair(nstride),idum,pi)

C (Original Comment) - Sort by z coordinates
c     zsort clusters 
            
            call z_heapsort(natoms,pair,qp)				! Cluster 1
            call z_heapsort(natoms,pair(nstride),qp(natoms+1))		! Cluster 2
            
C (Original Comments)
c     ** mates selected clusters **
c     
c     mate_1pt_random - one point crossover about a random position
c     mate_1pt_weighted - one point crossover about a random position
c     mate_2pt - two point crossover about two randomly chosen positions
         
C            write(*,'('' DEBUG '',/)')
  
            if ( nels .eq. 2 ) then
C If bimetallic

               if (mtype .eq. 1) then
C Anyway, random mating
                  call mate_1pt_multi_random(natoms,pfunc,
     &                 pop((nclust+j-1)*nstride+1),
     &                 q((nclust+j-1)*natoms+1),pair,qp,idum)
               else if (mtype .eq. 2) then
C Or weighted crossover
                  call mate_1pt_multi_weighted(natoms,pfunc,
     &                 pop((nclust+j-1)*nstride+1),
     &                 q((nclust+j-1)*natoms+1),pair,qp,
     &                 fit(str(1)),fit(str(2)))
               else
C 2pt works fine. This is not default though
                  call mate_2pt_multi_random(natoms,pfunc,
     &                 na,nb,pop((nclust+j-1)*nstride+1),
     &                 q((nclust+j-1)*natoms+1),pair,qp,
     &                 idum)
               end if
            else
C (Original Comments)
c     one and two point crossover for single element clusters
               if (mtype .eq. 1) then
C 1 point random
                  call mate_1pt_random(natoms,
     &                 pop((nclust+j-1)*nstride+1),pair,idum)
               else if (mtype .eq. 2) then
C 1 point weighted
                  call mate_1pt_weighted(natoms,
     &                 pop((nclust+j-1)*nstride+1),pair,fit(str(1)),
     &                 fit(str(2)))
               else
C Two point 
                  call mate_2pt(natoms,pop((nclust+j-1)*nstride+1),pair,
     &                 idum)
               end if
            end if
       
C            write(*,'('' DEBUG '',/)')      
      
CC MUTATION CC
C (Original Comments)
c       increment counter for mutant array
            if ( write_stats ) mutant(nclust+j) = 1	! Edit mutant statistics
            
            if (mscheme .eq. 1) then
C If mutation is in place               
               mfact=ran3(idum)		! Random Number
               
               if (mfact .lt. mrate) then 
C If factor is greater than mutation rate
C Debug
C                  write(*,'(''Mutating'')')

C (Original Comments)                  
c     ** mutation **
c     mutation scheme 1 - mutate offspring only (and not a copy!)
C All mutantion functions are in mutate.f

C If mutate_exchange defined, do this
                  if (muttype .eq. 4) then
C Mutate Exchange
                     call mutate_exchange(natoms,na,nb,
     &                    pop((nclust+j-1)*nstride+1),
     &                    q((nclust+j-1)*natoms+1),
     &                    idum,pfunc,meswaps)
                  else
C Otherwise, check if a mutate_Exchange rate is defined
C If so pull random number and check probability of mutation types
                     if (merate .gt. 0.d0) then
                        mfact=ran3(idum)         ! Random Number

                        if (mfact .lt. merate) then
C Debug
C                          write(*,'(''Swapping Atoms'')')
C Probability is within region, do mutate exchange
                          call mutate_exchange(natoms,na,nb,
     &                    pop((nclust+j-1)*nstride+1),
     &                    q((nclust+j-1)*natoms+1),
     &                    idum,pfunc,meswaps)
                        else
C Debug
C                          write(*,'(''Moving Atoms'')')
C Otherwise do other type of mutation as this is more probable
                          if (muttype .eq. 1) then
C Mutate replace
                            call mutate_replace(natoms,na,nb,qa,qb,r_ab,
     &                      pop((nclust+j-1)*nstride+1),
     &                      q((nclust+j-1)*natoms+1),idum,range,nels)
                          else if (muttype .eq. 2) then
C Mutate rotate
                  call mutate_rotate(natoms,pop((nclust+j-1)*nstride+1),
     &                      q((nclust+j-1)*natoms+1),idum,pi)
                          else if(muttype .eq. 3) then
C Mutate Move
                   call mutate_move(natoms,pop((nclust+j-1)*nstride+1),
     &                      idum,range)
                          end if
C End of possible link with mutate_Exchange
                        end if
                     else
C Go through different options and do which is chosen
                       if (muttype .eq. 1) then
C Mutate replace
                         call mutate_replace(natoms,na,nb,qa,qb,r_ab,
     &                    pop((nclust+j-1)*nstride+1),
     &                    q((nclust+j-1)*natoms+1),idum,range,nels)
                       else if (muttype .eq. 2) then
C Mutate rotate
                  call mutate_rotate(natoms,pop((nclust+j-1)*nstride+1),
     &                    q((nclust+j-1)*natoms+1),idum,pi)
                       else if(muttype .eq. 3) then
C Mutate Move
                   call mutate_move(natoms,pop((nclust+j-1)*nstride+1),
     &                    idum,range)
                       end if
C All mutation options complete
                     end if
                  end if
                  
                  if ( write_stats ) mutant(nclust+j) = 2

               end if
       
c     relax offspring
      
               call minimiser(natoms,na,nb,idum,pfunc,idamp,ihard,nels,
     &             pop((nclust+j-1)*nstride+1),q((nclust+j-1)*natoms+1),
     &             alpha,r_ab,qa,qb,B_ab,rho_ab,B_bb,rho_bb,de,re,a2,a3,
     &             coeff,rcut2,range,nspec,aij,vij,pij,qij,r0,
     &             c1_sisi,c2_sisi,c3_sisi,c4_sisi,
     &             c1_sio,c2_sio,c3_sio,c4_sio,
     &             c1_oo,c2_oo,c3_oo,c4_oo,
     &             ru_sisi,ru_sio,ru_oo,cutoff_gupta,cutoff_start,
     &             cutoff_end,a3a,a4a,a5a,x3x,x4x,x5x,i,j,
     &             namea,nameb,mass_a,mass_b)
   
C Save energy values
               energy(nclust+j) = pop((nclust+j)*nstride)
C Increment number of minimised structures
               nmin = nmin + 1
C Check if we have hit target energy
               if ( dabs(energy(nclust+j)-gmin) .lt. 2.d-6 
     &              .and. .not. gflag) then
                  write(*,'(a,1x,i4)')
     &                 'GMIN found by minimization number ',nmin
                  gflag = .true.
               end if 
C (Original Comment)
c     required for energy sort and other operators
               nmut = 0
               
            end if
         end do

C         write(*,'('' DEBUG '',/)')

C (Original Comments)
c     mutation scheme 2 - mutate copies (of population and offspring)

         if (mscheme .eq. 2) then
             nmut = 0
            
            do j=1,nclust+noff
C For all clusters
               
               mfact=ran3(idum)		! Random number
               
               if (mfact .lt. mrate) then 
C Debug
C                  write(*,'(''Mutating'')')
C If we meet the mutation criteria
                  
c     ** mutation **
c     All mutation schemes stored in mutate.f     
                  
                  do k=1,3*natoms
C For all coordinates, move values back in population to make space
                     pop((nclust+noff+nmut)*nstride+k) =
     &                    pop((j-1)*nstride+k)
                     ! print *,pop((nclust+noff+nmut)*nstride+k)
                  end do
                  
                  do k=1,natoms
C For all clusters, move charges back to make space for new clusters
                     q((nclust+noff+nmut)*natoms+k) = 
     &                    q((j-1)*natoms+k)
                     ! print *,q((nclust+noff+nmut)*natoms+k)
                  end do

C If mutate_exchange defined, do this
                  if (muttype .eq. 4) then
C Mutate Exchange
                     call mutate_exchange(natoms,na,nb,
     &               pop((nclust+noff+nmut)*nstride+1),
     &               q((nclust+noff+nmut)*natoms+1),
     &               idum,pfunc,meswaps)
                  else
C Otherwise, check if a mutate_Exchange rate is defined
C If so pull random number and check probability of mutation types
                     if (merate .gt. 0.d0) then
                        mfact=ran3(idum)         ! Random Number

                        if (mfact .lt. merate) then
C Debug
C                          write(*,'(''Swapping Atoms'')')
C Probability is within region, do mutate exchange
                          call mutate_exchange(natoms,na,nb,
     &                    pop((nclust+noff+nmut)*nstride+1),
     &                    q((nclust+noff+nmut)*natoms+1),
     &                    idum,pfunc,meswaps)
                        else
C Debug
C                          write(*,'(''Moving Atoms'')')
C Otherwise do other type of mutation as this is more probable
                          if (muttype .eq. 1) then
C Mutate replace
                            call mutate_replace(natoms,na,nb,qa,qb,r_ab,
     &                      pop((nclust+noff+nmut)*nstride+1),
     &                   q((nclust+noff+nmut)*natoms+1),idum,range,nels)
                          else if (muttype .eq. 2) then
C Mutate rotate
                            call mutate_rotate(natoms,
     &                      pop((nclust+noff+nmut)*nstride+1),
     &                      q((nclust+noff+nmut)*natoms+1),idum,pi)
                          else if(muttype .eq. 3) then
C Mutate Move
                            call mutate_move(natoms,
     &                     pop((nclust+noff+nmut)*nstride+1),idum,range)
                          end if
C End of possible link with mutate_Exchange
                        end if
                     else
C Go through different options and do which is chosen
                       if (muttype .eq. 1) then
C Mutate replace
                         call mutate_replace(natoms,na,nb,qa,qb,r_ab,
     &                   pop((nclust+noff+nmut)*nstride+1),
     &                   q((nclust+noff+nmut)*natoms+1),idum,range,nels)
                       else if (muttype .eq. 2) then
C Mutate rotate
                         call mutate_rotate(natoms,
     &                   pop((nclust+noff+nmut)*nstride+1),
     &                   q((nclust+noff+nmut)*natoms+1),idum,pi)
                       else if(muttype .eq. 3) then
C Mutate Move
                          call mutate_move(natoms,
     &                    pop((nclust+noff+nmut)*nstride+1),idum,range)
                       end if
C All mutation options complete
                     end if
                  end if

CCCCCCCCCCCCCCCCCCCCCCCC
C
C                 write(*,'(''Mutated'')')
C
C                  do k=1,3*natoms
C                     print *,pop((nclust+noff+nmut)*nstride+k)
C                  end do
C
C                  do k=1,natoms
C                     print *,q((nclust+noff+nmut)*natoms+k)
C                  end do
C
CCCCCCCCCCCCCCCCCCCCCCCC

C Count mutants                  
                  nmut = nmut + 1	! Number of Mutants ++

C                  print *,nmut
                  
                  if (write_stats) mutant(nclust+noff+nmut) = 2
 
               end if
            end do
            
            do j=nclust+1,nclust+noff+nmut
C For all new clusters relax offspring

C      write(*,*)i      ! Debug

C And minimise structures               
               call minimiser(natoms,na,nb,idum,pfunc,idamp,ihard,nels,
     &            pop((j-1)*nstride+1),q((j-1)*natoms+1),
     &            alpha,r_ab,qa,qb,B_ab,rho_ab,B_bb,rho_bb,de,re,a2,a3,
     &            coeff,rcut2,range,nspec,aij,vij,pij,qij,r0,
     &            c1_sisi,c2_sisi,c3_sisi,c4_sisi,
     &            c1_sio,c2_sio,c3_sio,c4_sio,
     &            c1_oo,c2_oo,c3_oo,c4_oo,
     &            ru_sisi,ru_sio,ru_oo,cutoff_gupta,cutoff_start,
     &            cutoff_end,a3a,a4a,a5a,x3x,x4x,x5x,i,j,
     &            namea,nameb,mass_a,mass_b)



c      write(*,*) 'energy before moisort',pop(j*nstride)


C Save energy               
               energy(j) = pop(j*nstride)
C Increment minimised structures               
               nmin = nmin + 1
C Check if we have the minimum we are looking for
               if ( dabs(energy(j)-gmin) .lt. 2.d-6 
     &              .and. .not. gflag) then
                  write(*,'(a,1x,i4)')
     &                 'GMIN found by minimization number ',nmin
                  gflag = .true.
               end if
               
            end do
         end if

C Moment of inertia calculation

          call moi(mass_a, mass_b, natoms, nstride, na, nb,
     &        pop, q, nclust+noff+nmut,
     &        inertia)

C Comparing moment of inertia for initial
C cluster to offspring/mutants


      call moisort(inertia, energy,nclust, noff, nmut, temp)
c         write(*,*) 'energyafter', energy 
 






C         write(*,'('' DEBUG '',/)')
      
C (Original Comments)
c       sort clusters by energy (best first)
C Subroutine in pop.f
   
         call esort(natoms,nstride,pop,q,energy,mutant,
     &        nclust+noff+nmut,nels)

C (Original Comments)
c       removes similar clusters
C Subroutine in pop.f
   
         if (remove) call maintain(natoms,nclust,nstride,nels,
     &        pop,q,energy,nclust+noff+nmut,pfunc,GAconv)

C (original Comments)
c     adds high energy mutant to population 
C Subroutine in pop.f
         
         if (highemut) call highE_mutant(natoms,nclust,nstride,
     &        noff+nmut,idum,pfunc,pop,q,mutant,energy)

C (Original Comments)
c     predator
C Subroutine in pop.f
         if (pred) call e_pred(natoms,nclust,nstride,noff+nmut,nels,pop,
     &     q,energy,epred)

C (Original Comments)
c     write out best/worst/average energy to screen
c     & summary file
C Subroutine in print_sum.f
         call print_sum(nclust,energy,i)

         if ( write_clust ) then
C (Original Comments)
c     write geometries of current generation to file
C Subroutine in print_sum.f
            call write_gen(natoms,nclust,nstride,pfunc,nels,pop,energy,
     &           q,el,namea,nameb)
         end if

C (Original Comment)
c     write out energies
C Subroutine in print_sum.f
         if ( write_energy) call write_en(nclust,i,energy)

C (Original Comment)
c     write out stats
         if ( write_stats ) then
            
            mut = 0
            mat = 0
            
            do j=1,nclust
C For all clusters
               if (mutant(j) .eq. 1) mat = mat + 1	! Check if mated
               if (mutant(j) .eq. 2) mut = mut + 1	! Check if mutated
            end do
            
C Output resulting stats
            write(*,'('' Number of offspring accepted = '',i2)')mat
          write(*,'('' Number of mutated offspring accepted = '',i2,/)')
     &           mut

          write(9,'(i3,2x,i2,2x,i2)')i,mat,mut

C Write out nearest neighbours as well (Sven's work)
C          call nextNeigh(natoms,nclust,nstride,nels,pop,q,
C     &                nclust+nmut+noff,r_ab,i)

       end if

c     end run of program if energy = target (given in gmin)

         if ( gmin .lt. 0 ) then
C Check if minimum energy criteria
            if ( dabs(energy(1)-gmin) .lt. 2.d-6 ) then
C See if cluster energy is below target

C Write coordinates to file
               call xyz(natoms,pfunc,pop,q,bscale,'end ',i,flen,
     &              namea,nameb,el,fname,nels)
               write(*,
     &'(''Target energy, '',f11.6,'', found during generation '',i3,/)')
     &              energy(1),i
               stop		! End program
            end if
         end if
         
C (Original Comments)
c     test for convergence - end if no change in population over tcount 
c     generations
         
         tflag = .false.
         if (pfunc .ne. 6) GAconv = 1.d-6  ! if not a DFT calculation
         call term(nclust,iterm,tcount,energy,tflag,GAconv,oldCluster)
         if (tflag) then
            write(*,'(''Convergence on generation: '',i3,/)')i
            write(*,'('' Writing xyz file ...... '',/)')

C Write output coordinates to file
            call xyz(natoms,pfunc,pop,q,bscale,
     &           'term ',i,flen,namea,nameb,el,fname,nels)
            stop		! End program
         end if
         
c     end of one generation    
         
      end do

C No convergence
C Write output coordinates to file 
      call xyz(natoms,pfunc,pop(1),q(1),bscale,
     &     'none ',0,flen,namea,nameb,el,fname,nels)

      stop
      
      end

