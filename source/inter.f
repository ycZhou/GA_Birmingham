C    Subroutine interpret:
C    This will read an input file, and assign variables
C    Commenting added 20/10/2010, Andrew Logsdail
C    as part of the on-going tidy up project!

      subroutine interpret(potential)

      use commons
      use potentials
      implicit none
C     Boolean values

C     String values
      character*80 line
      character*20 potential 	! Potential file name

C     String variables with lengths defined in ga.param
      character*(tokbufflen) buffer	! String buffer
      character*(maxtoklen) token 	! Holds our token value

C     Internal int variables
      integer i
      integer j
      integer k
      integer last
      integer n
      integer tok			! Current token
      integer tlen
      integer ttok			! Total number of tokens
      integer chartoint			! Integer representation of character
      integer lent			! Length of string
      integer ln
      integer number			! Input number

C     Integer arrays
      integer lno(maxtokens)
      integer bloc(maxtokens)
      integer blen(maxtokens)
      integer loc(maxtokens)
      integer len(maxtoklen)

C     Internal routine character variables
      real*8 chartodb
      real*8 dnum

C     Boolean routine variables
C     These are all used to see if the variables have been set
      logical lprocessors		!
      logical lnatoms			!
      logical lnatoms_a			!
      logical lnatoms_b			!
      logical lnclust			!
      logical lngen			!
      logical lnoff			!
      logical lmrate			!
      logical lgmin			!
      logical lterm			!
      logical ltsize			!
      logical lptype			!
      logical lmate			!
      logical lsel			!
      logical lfit			!
      logical lmut			!
      logical lfile			!
      logical lseed			!
      logical error			!
      logical lmscheme			!
      logical lels			!
      logical lmeswaps			!
      logical lmerate			!

C     External functions
      external chartoint		! Character to int
      external chartodb			! Character to double

C     (Original comment)
c     read and tokenise input file

      call tokeniser(ttok,tlen,blen,bloc,lno,buffer)
      
C     (Original comment)
c     set all flags to false
C     If set at the end we get a print out using these variables
C     Firstly all the internal (subroutine booleans)

      lprocessors = .false.
      lnatoms   = .false.
      lnatoms_a = .false.
      lnatoms_b = .false.
      lels      = .false.
      lnclust   = .false.
      lngen     = .false.
      lnoff     = .false.
      lmrate    = .false.
      lgmin     = .false.
      lterm     = .false.
      ltsize    = .false.
      lptype    = .false.
      lmate     = .false.
      lsel      = .false.
      lfit      = .false.
      lmut      = .false.
      lfile     = .false.
      lseed     = .false.
      lmeswaps  = .false.
      lmerate   = .false.

C     These are program variables
      highemut  = .false.
      remove    = .false.
      write_clust = .false.
      write_energy = .false.
      write_stats = .false.
      pred = .false.
      error = .false.
      lmscheme = .false.
C     Added by Andrew Logsdail 2009
      cutoff_gupta = .false.
      restart = .false.

C     (Original comment)
c     loop until all tokens processed

C     Set tok at starting point
      tok = 1
C     While less than total tokens...
      do while ( tok .le. ttok )	! This will make sure we iterate over the whole file

C Debug
C         write(*,*)tok

C        Let's clear the token value each time to make sure we get true values
         token = ''
C        Get next token, assign to "token" with length "lent"
         call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok) 	! Get file
C         print *,token(1:lent)						! Print to screen (check)

C     (Original comments)          
c     begin large number of ifs
c     start with parameters expecting integer numbers

C        Number of atoms (monometallic)
         if ( token .eq. 'natoms' ) then
            
C           Expect int will read the next token as an integer, and assign to number
C           Works consistently over all integer values at present
            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lnatoms,error,ln)		! test() is used to check the value is allowed. This is used throughout this code

            lnatoms = .true.
            natoms = number

C        Number of atoms (type a)
         else if ( token .eq. 'natoms_a' ) then
            
            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lnatoms,error,ln)
            
            lnatoms_a = .true.
            na = number
            
C        Number of atoms (type b)
         else if ( token .eq. 'natoms_b' ) then
            
            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lnatoms,error,ln)

            lnatoms_b = .true.
            nb = number

C        Number of clusters in population
         else if ( token .eq. 'nclust') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lnclust,error,ln)
            
            lnclust = .true.
            nclust = number

C        Number of generations
         else if (token .eq. 'ngen') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lngen,error,ln)

            lngen = .true.
            ngen = number
            
C        Number of offspring
         else if (token .eq. 'noff') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lnoff,error,ln)

            lnoff = .true.
            noff = number

C        Terminating term
         else if (token .eq. 'term') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lterm,error,ln)

            lterm = .true.
            tcount = number

C         Tournament size for mating selections
          else if (token .eq. 'tsize') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(ltsize,error,ln)

            ltsize = .true.
            tsize = number  

C        Check for starting seed value for random numbers
         else if (token .eq. 'seed') then
            
            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)          
            call test(lseed,error,ln)

            lseed = .true.
            idum = number

C        Number of elements (mono- or bimetallic)
         else if (token .eq. 'nels') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lels,error,ln)

            lels = .true.
            nels = number

C	Get number of processors for OMP
         else if (token .eq. 'processors') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,number,
     &           error)
            call test(lprocessors,error,ln)

            lprocessors = .true.
            processors = number

C     (Original comment)
c     move on to values expecting real numbers
C     OK - so here we are going to take in decimal numbers

C        Expect mutation rate as percent of total population
         else if (token .eq. 'mrate') then

C           Read next token as double
            call expect_dble(token,ttok,tok,buffer,blen,bloc,lno,dnum,
     &           error)
            call test(lmrate,error,ln)
            
            lmrate = .true.
            mrate = dnum

C        Energy minimum
C        Program will terminate if cluster is found with energy below this threshold
         else if (token .eq. 'gmin') then

            call expect_dble(token,ttok,tok,buffer,blen,bloc,lno,dnum,
     &           error)
            call test(lgmin,error,ln)
            
            lgmin = .true.
            gmin = dnum

C        High energy predator as percentage of total population
         else if (token .eq. 'predator') then

            call expect_dble(token,ttok,tok,buffer,blen,bloc,lno,dnum,
     &           error)
            call test(pred,error,ln)
            
            pred = .true.
            epred = dnum

C        Now all the numerics have been read.
C        Next up we will read in the different poptentials

C        Morse potential
         else if (token .eq. 'morse') then
               
            call test(lptype,error,ln)

            tok = tok + 1
            lptype = .true.
            ptype = 1

C        All the Murrell Mottram's should be compiled together into one if statement,
C        similar to the gupta potential. This will make the program more dynamic.
C        Perhaps something for a rainy day...
C        ..... or sunny day. Same day completion

C        Here we just identify that we are dealing with Murrell Mottram, and save the
C        token as "potential".  
C        !!!!!!!!!!!!!!!!!!!!!!!
C        Separating out all the variables appears to have solved this problem
C        Think it was just a memory issue with things being rewritten
C        !!!!!!!!!!!!!!!!!!!!!!!
        else if (token(1:3) .eq. 'mm_') then
C         Unlike the others we need to specify the substring here

C            Error test          
C            print *,fname

            call test(lptype,error,ln)

C            Error test
C            print *,fname

            potential = token
            tok = tok + 1
            lptype = .true.
            ptype = 2                  ! All MM type potentials are type 2

C        Magnesium Oxide Potential
C        Though the form is not obvious from here, perhaps it will show later in the code
         else if (token .eq. 'MgO') then
            
            call test(lptype,error,ln)
            
            tok = tok + 1
            lptype = .true.
            ptype = 7			! Potential type 7 is MgO
            
C        Zinc Oxide Potential
C        Same as above, form is not obvious
         else if (token .eq. 'ZnO') then
            
            call test(lptype,error,ln)
            
            tok = tok + 1
            lptype = .true.
            ptype = 8			! Potential type 8 is ZnO

C        ?? Unknown potential - generic ion ??
C        Again form not obvious. This needs more research!
         else if (token .eq. 'gen_ion') then
            
            call test(lptype,error,ln)
            
            tok = tok + 1
            lptype = .true.
            ptype = 9			! Potential type 9
            
C        Gupta potential - deals with mono- and bi-metallic
C        Here we just identify that we are dealing with gupta, and save the
C        token as "potential".  
C        !!!!!!!!!!!!!!!!!!!!!!!
C        In gfortran this command seems to wipe the fname variable.
C        Why??
C        !!!!!!!!!!!!!!!!!!!!!!!
C        Separating out all the variables appears to have solved this problem
C        Think it was just a memory issue with things being rewritten
C        !!!!!!!!!!!!!!!!!!!!!!!
        else if (token(1:6) .eq. 'gupta_') then
C         Unlike the others we need to specify the substring here
 
C            Error test          
C            print *,fname
 
            call test(lptype,error,ln)

C            Error test
C            print *,fname

            potential = token
            tok = tok + 1
            lptype = .true.
            ptype = 10			! All gupta type potentials are type 10

C        Silicon Oxide
C        Form not obvious
         else if (token .eq. 'SiO_TTAM') then

            call test(lptype,error,ln)

            tok = tok + 1
            lptype = .true.
            ptype = 24			! Potential type is 24

C        Silicon Oxide
C        Form not obvious
         else if (token .eq. 'SiO_FB') then

            call test(lptype,error,ln)

            tok = tok + 1
            lptype = .true.
            ptype = 25			! Potential type 25

C        DFT Calculation 
         else if (token(1:4) .eq. 'DFT_') then

            call test(lptype,error,ln)

            potential = token
            tok = tok + 1
            lptype = .true.
            ptype = 26                  ! Potential type 26 for all DFT

C        That is all the potential types read
C     (Original comment)
c     mating

C        Read in mating type for populations
C        Firstly 1pt mating, random
         else if (token .eq. 'mate_1pt_random') then

            call test(lmate,error,ln)
            
            tok = tok + 1
            mtype = 1			! Mating type 1
            lmate = .true.
            
C        1pt mating, weighted
         else if (token .eq. 'mate_1pt_weighted') then
            
            call test(lmate,error,ln)
            
            tok = tok + 1
            mtype = 2			! Mating type 2 (Default)
            lmate = .true.
            
C        2pt mating
         else if (token .eq. 'mate_2pt') then
            
            call test(lmate,error,ln)
            
            tok = tok + 1
            mtype = 3			! Mating type 3
            lmate = .true.

C        1pt mating, random, separated
         else if (token .eq. 'mate_1pt_random_sep') then

            call test(lmate,error,ln)
            
            tok = tok + 1
            mtype = 4			! Mating type 4
            lmate = .true.

C        ?? Unknown mating type ??
C        This needs further investigation to be commented well
         else if (token .eq. 'mate_inter') then

            call test(lmate,error,ln)
            
            tok = tok + 1
            mtype = 5			! Mating type 5
            lmate = .true.
            
C     (Original Comment)
c     selection

C        Roulette selection for mating
         else if (token .eq. 'roulette') then
            
            call test(lsel,error,ln)
            
            tok = tok + 1
            stype = 1			! Selection type 1
            lsel = .true.

C        Tournament selection for mating
         else if (token .eq. 'tournament') then
            
            call test(lsel,error,ln)
            
            tok = tok + 1
            stype = 2
            lsel = .true.
            
C     (Original comment)
c     fitness

C        Linear fitness function
         else if (token .eq. 'linear_fitness') then
            
            call test(lfit,error,ln)
            
            tok = tok + 1
            ftype = 1			! Fitness type 1
            lfit = .true.
            
C        Power fitness
         else if (token .eq. 'power_fitness') then
            
            call test(lfit,error,ln)
            
            tok = tok + 1
            ftype = 2			! Fitness type 2
            lfit = .true.

C	 Exponetial fitness
         else if (token .eq. 'exp_fitness') then
            
            call test(lfit,error,ln)
            
            tok = tok + 1            
            ftype = 3			! Fitness type 3
            lfit = .true.
            
C        Tanh fitness
         else if (token .eq. 'tanh_fitness') then
            
            call test(lfit,error,ln)
            
            tok = tok + 1
            ftype = 4			! Fitness type 4 (Default)
            lfit = .true.

C        Cos2 fitness function
         else if (token(1:lent) .eq. 'cos2_fitness') then
            
            call test(lfit,error,ln)
            
            tok = tok + 1
            ftype = 5			! Fitness type 5
            lfit = .true.

C     (Original comment)
c     mutation

C        This expects an integer value
C        Presumably we can specify mutation scheme directly by number?
C        Needs testing from the command line.
         else if(token .eq. 'mutation_scheme') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,
     &        number,error)
            call test(lmscheme,error,ln)

            lmscheme = .true.
            mscheme = number

C        Mutate replace to be used
         else if (token .eq. 'mutate_replace') then
            
            call test(lmut,error,ln)
            
            tok = tok + 1
            muttype = 1			! Mutation Type 1
            lmut = .true.

C        Use mutate and rotate (default)
         else if (token .eq. 'mutate_rotate') then
            
            call test(lmut,error,ln)
            
            tok = tok + 1
            muttype = 2			! Mutation Type 2
            lmut = .true.

C        Use mutate and move
         else if (token .eq. 'mutate_move') then 
            
            call test(lmut,error,ln)
            
            tok = tok + 1
            muttype = 3			! Mutation Type 3
            lmut = .true.
            
C        Exchange atoms as part of mutation process
         else if (token .eq. 'mutate_exchange') then
            
            call test(lmut,error,ln)

            tok = tok + 1
            muttype = 4			! Mutation Type 4
            lmut = .true.

CC Added November 2010
C        This expects an integer value
C        And defines the number of bimetallic position swaps we'd use
C        in accompaniment with other mutation schemes
         else if(token .eq. 'me_swaps') then

            call expect_int(token,ttok,tok,buffer,blen,bloc,lno,
     &        number,error)
            call test(lmeswaps,error,ln)

            lmeswaps = .true.
            meswaps = number

C        Here we define the rate at which mutate exchange happens ON TOP of any other mutation technique
C        For everytime, set to 1, for never set to 0
C        Remember that when mutate exchange is explicitly specified this will have no influence
         else if (token .eq. 'me_rate') then

C           Read next token as double
            call expect_dble(token,ttok,tok,buffer,blen,bloc,lno,dnum,
     &           error)
            call test(lmerate,error,ln)

            lmerate = .true.
            merate = dnum
            
C     (Original comments)
c     high energy mutants
            
C        Allow high energy mutants
         else if (token .eq. 'high_energy_mutants') then
            
            tok = tok + 1 
            highemut = .true.
            
C        remove clusters of similar energy
         else if (token .eq. 'remove_similar' 
     &                  .or. token .eq. 'remove') then
            
            tok = tok + 1
            remove = .true.

C        write out all clusters
         else if (token .eq. 'write_clusters') then

            tok = tok + 1
            write_clust = .true.

C        write out all energies
         else if (token .eq. 'write_energies') then

            tok = tok + 1
            write_energy = .true.

C     Andrew Logsdail additions to the code, 2009
C     Here I am introduccing the application of the gupta_cutoff if wanted.
C     Command name will be called "cutoff_gupta"

         else if (token .eq. 'cutoff_gupta') then
            tok = tok + 1
            cutoff_gupta = .true.

C    This will need to be accompanied by cutoff_parameter file, as is used for the Minimiser
C    This will contain cutoff_start, cutoff_end, a3, a4, a5, x3, x4 and x5
C    End Andrew Logsdail additions

C Adding functionality to do restarts, 2011
         else if (token .eq. 'restart') then
            tok = tok + 1
            restart = .true.

C     write out number of offspring and mutants accepted 
         else if (token .eq. 'write_stats') then

            tok = tok + 1
            write_stats = .true.

C     (Original comment)
c     file name

C        This gives us the preceding text to all filenames
         else if (token .eq. 'file') then
                  
C           Get next token, check if it is an equals sign
            tok = tok + 1			! Move to next token
            call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,
     &              tok)

            if (token(1:lent) .eq. '=') then

C              Get next token, check it is not too long.
               call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,
     &              tok+1)

               if (lent .ge. maxtoklen) then
                  write(*,
     & '(''ERROR on line '',i2,'' filename too long'')')ln
                  error = .true.
               end if

C              !! Fixed problem !! Here in case of future requirement
C              Tests to figure out why we are losing the filename
C               print *,token
C               print *,fname

C              Ok, now we transfer the fname.
C              Why is this done character by character not as transfer of string?

C               do i=1,lent
C                  fname(i:i) = token(i:i)
C               end do

C              Let's try just as a direct assignment, more updated code...
               fname = token
C              Works fine. Think the previous code was to prevent cross-contamination
C              but we've negated this by cleaning the token variable each time in the loop

C              !! FIXED PROBLEM !! Here in case of future requirements
C              Test to figure out why we are losing the filename
C               print *,token
C               print *,fname
               
               flen = lent			! Remember length of filename
               tok = tok + 2			! Move along three places as we've pulled extra bits of information
               lfile = .true.			! Loaded file type is true

            else

C              Print error if filename is not found
               write(*,
     & '(''ERROR on line '',i2,'' filename expected'')')ln
               tok = tok + 1			! Move along 2 places for "file" and "="

            end if

C        Final else statement in this collection...
         else

C            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C            unrecognised token, throw ERROR
C            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             write(*,
     &     '(''ERROR on line '',i2,'', unknown parameter '')')ln
            error = .true.
               
            tok = tok + 1

         endif

      end do
C
C     END OF INPUT FILE INTERPRETATION
C

C     Check parameters. Have we got...
C     A potential?
      if (.not. lptype) then
         write(*,
     & '(''ERROR, type of potential not specified'')')
         error = .true.
      end if

c     Warning: gupta_cutoff useless if not used for gupta potential
      if (ptype .eq. 10. .and. cutoff_gupta) then
            write(*,88)
 88         format('Warning, Cutoff not used with any other'
     &           ,' potential than the Gupta')
      end if


C     A seed for random numbers?
      if (.not. lseed) then
         write(*,
     & '(''ERROR, random number seed not specified'')')
         error = .true.
      end if

C     Number of elements value?
      if (.not. lels) then
         write(*,
     & '(''ERROR, Number of different atom types not specified'')')
         error = .true.
      end if

C     If monometallic are number of atoms less than max and greater than 2 (i.e. a dimer)
      if (nels .eq. 1) then

         if (lnatoms) then
            if (natoms .gt. nmax .or. natoms .lt. 2) then
               write(*,
     &     '(''ERROR, number of atoms must be between 2 and '',i3)')
     &              nmax
               error = .true.
            end if
C        Else atoms are not specificed
         else
            write(*,
     &           '(''ERROR, number of atoms not specified'')')
            error = .true.
         end if

      end if

c     Check if number of clusters is defined?
      if (lnclust) then
         if (.not. lnoff) noff = 0.8*nclust		! If number of offspring is undefined then take it as 80%
         if (nclust+noff .gt. ncmax) then		! Check total is not greater than allowed total
            write(*,
     & '(''ERROR, number of clusters is too large, max = '',i3)')
     &           ncmax
            write(*,
     & '(4x,''This is the total no of clusters: nclust + noff'')')
            if (.not. lnoff) then			! Inform user default offspring setting used
               write(*,
     & '(''Used default value for noff of 0.8*nclust, = '',i2)')noff
            end if
            error = .true.
         end if
         if (nclust .le. 2) then			! Check number of clusters constitutes a population
            write(*,
     & '(''ERROR, number of clusters must be greater than 2'')')
         end if
      else
         write(*,
     & '(''ERROR, number of clusters not specified'')')
         error = .true.
      end if

c     Check number of generations is specficed?
      if (.not. lngen) then
         write(*,
     & '(''ERROR, number of generations not specified'')')
         error = .true.
      end if

c     test tournament selection variables...
      if (lsel .and. stype .eq. 2) then							! If selection type defined and equals = 2
         if (.not. ltsize) then							! If size of tournament is not defined
            write(*,
     & '(''ERROR, no tournament size given for tournament selection'')')
            error = .true.
         end if
         if (ltsize .and. tsize .gt. nclust) then					! Check tournament size isn't greater than number of clusters
             write(*,
     & '(''ERROR, tournament size larger than population size'')')
             error = .true.
         end if
      else
         if (ltsize) then								! Arbitrary warning for roulette selection
            write(*,4)
 4          format('Warning, specifying a tournament size when'
     &              ,' using roulette selection is meaningless')
         end if
      end if

c     Warning: selection of fitness type if using tournament selection
      if (stype .eq. 2. .and. lfit) then
            write(*,33)
 33         format('Warning, fitness relationship not used with'
     &           ,' tournament selection')
      end if

c     Check mutation rate...
      if (lmrate) then
         if (mrate .gt. 1.d0 .or. mrate .lt. 0.d0) then		! Check mutation rate fits is percentage
            write(*,
     & '(''ERROR, mutation rate should lie in the range 0.0-1.0'')')
            error = .true.
         end if
      end if

c     Check mutation rate for mutate_exchange ...
      if (lmerate) then
         if (merate .gt. 1.d0 .or. merate .lt. 0.d0) then         ! Check mutation exchange rate fits is percentage
            write(*,
     & '(''ERROR, mutation rate should lie in the range 0.0-1.0 
     & for mutate_exhange probability'')')
            error = .true.
         else if (muttype .eq. 4) then
            write(*,'(''Mutation rate probability for mutate
     &  exchange coupled with mutation_exchange serves no function.
     &  Continuing anyway...'')')
         end if
      end if


c     Mutation scheme - choices are only 1 and 2
      if (lmscheme) then
         if (mscheme .ne. 1 .and. mscheme .ne. 2) then
            write(*,
     & '(''ERROR, invalid mutation scheme; choices are 1 and 2'')')
            error = .true.
         end if
      end if

c     Check number of atoms of each type for multi element clusters
      if (nels .eq. 2) then

C        Check number of type a is defined
         if (.not. lnatoms_a) then
            write(*,
     &  '(''ERROR, number of atoms of type a not specified'')')
            error = .true.
         end if
C        Check number of type b
         if (.not. lnatoms_b) then
            write(*,
     & '(''ERROR, number of atoms of type b not specified'')')
            error = .true.
         end if
         if (lnatoms_a .and. lnatoms_b) then
C           Check that total number of atoms is within limits
c           Sum na and nb to give total atoms
            natoms = na + nb
C           Then compare
            if (natoms .gt. nmax .or. natoms .lt. 2) then
               write(*,
     &'(''ERROR, total number of atoms must lie in range 2 to '',i3)')
     &              nmax
               error = .true.
            end if
         end if

CC Redundant error
CC Mating type 2pt - does it work for bimetallics? We'll soon find out...
c     mate_2pt doesn't work for dual element clusters
C         if (lmate .and. mtype .eq. 3) then
C            write(*,
C     & '(''ERROR, no mate_2pt routine for dual element clusters'')')
C            error = .true.
C         end if

C     End checks for bimetallic system
      else
C     Checks for monometallic system

c     mutate_exchange meaningless for single element clusters
         
         if (lmut .and. muttype .eq. 4) then
            write(*,5)
 5             format('ERROR, mutate_exchange pointless for single 
     &              element clusters')
            error = .true.
         else if (lmeswaps .or. lmerate) then
            write(*,5)
            error = .true.
         end if

      end if

c     End of errors?
C     Let's print end of file if this a problem has occurred. 
      if (error) then
c     stop here if there are any errors

         write(*,'(/,''*******************************************'')')
         write(*,'(''   Run cancelled due to errors in input    '')')
         write(*,'(''*******************************************'',/)')

         stop
C     Otherwise print all information to screen
C     This is very verbose so I won't comment much of it
C     It should all be intuitive to read    
      else
   
 10      format('Number of atoms = ',i3)
 11      format('Number of atoms of type ',a,' = ',i3)

         write(*,'(''Input OK, continuing with GA run'',/)')

C     (Original comment)
c     start with input values
         
         if (ptype .eq. 1) then
            write(*,'(''Using Morse potential'')')
            write(*,10)natoms					! This prints number of atoms to screen, as per format 10
         else if (ptype .eq. 2) then
            write(*,'(''Using Murrell-Mottram potential'')')
            write(*,10)natoms
C         else if (ptype .eq. 3) then
C            write(*,'(''Using Murrell-Mottram Aluminium potential'')')
C            write(*,10)natoms
C         else if (ptype .eq. 4) then
C            write(*,'(''Using Murrell-Mottram Ytterbium potential'')')
C            write(*,10)natoms
C         else if (ptype .eq. 5) then
C            write(*,'(''Using Murrell-Mottram Lead potential'')')
C            write(*,10)natoms
C         else if (ptype .eq. 6) then
C          write(*,'(''Using user supplied Murrell-Mottram potential'')')
C            write(*,10)natoms
         else if (ptype .eq. 7) then
            write(*,'(''Using Rigid Ion potential for MgO'')')
            write(*,11)'a',na
            write(*,11)'b',nb
         else if (ptype .eq. 8) then
            write(*,'(''Using Rigid Ion potential for ZnO'')')
            write(*,11)'a',na
            write(*,11)'b',nb
         else if (ptype .eq. 9) then
            write(*,'(''Using user suppled Rigid Ion potential'')')
            write(*,11)'a',na
            write(*,11)'b',nb
         else if (ptype .eq. 10) then
            write(*,'(''Using Gupta Potential'')')
            write(*,11)'a',na
            write(*,11)'b',nb
         else if (ptype .eq. 24) then
            write(*,'(''Using TTAM potential for silica'')')
            write(*,11)'a',na
            write(*,11)'b',nb
         else if (ptype .eq. 25) then
            write(*,'(''Using FB potential for silica'')')
            write(*,11)'a',na
            write(*,11)'b',nb
         else if (ptype .eq. 26) then
            write(*,'(''Using DFT Calculations'')')
            write(*,11)'a',na
            write(*,11)'b',nb
C         else if (ptype .eq. 27) then
C            write(*,'(''Using DFT_RESTART option'')')
C            write(*,11)'a',na
C            write(*,11)'b',nb
         end if

C Added by Andrew Logsdail, 2009
C Noting of Gupta Cutoff Use
         if (cutoff_gupta) then
            write(*,'(''Using Gupta Cutoff '')')
         end if
C End of Logsdail additions
C Added by Andrew Logsdail, 2009
C Noting of restart
         if (restart) then
            write(*,'(''Restarting from previous coordinates '')')
         end if
C End of Logsdail additions


         write(*,'(''Number of clusters = '',i3)')nclust
         write(*,'(''Number of generations = '',i3)')ngen
         write(*,'(''Seed number = '',i3)')idum
        
C        Offspring 
         if (lnoff) then
            write(*,'(''Number of offspring = '',i3)')noff
         else
            write(*,12)noff
 12         format('Using DEFAULT value for number of offspring of'
     &           ,' 0.8*nclust = ',i2)
         end if

C        Mutation Scheme 
         if (lmscheme) then
            if (mscheme .eq. 1) then
               write(*, '("Using mutation scheme 1")')
            else if (mscheme .eq. 2) then
               write(*,'("Using mutation scheme 2")')
            end if
         else
C           Default mutation scheme
            mscheme = 2
            write(*,'("Using DEFAULT mutation scheme ",i1)')
     &           mscheme
         end if

C        Mutation Rate
         if (lmrate) then
            write(*,'(''Mutation rate = '',f4.2)')mrate
         else
c           Default mutation rate
            mrate = 0.1d0
            write(*,
     &   '(''Using DEFAULT value for mutation rate of '',f4.2)')mrate
         end if

C        Termination Counter
         if (lterm) then
            write(*,13)tcount
 13         format('GA will terminate if no change in population' 
     &           ,' over ',i3,' generations')
         else
c        Default termination counter
            tcount = ngen+1
            write(*,14)tcount
 14         format('Using DEFAULT termination counter of ',i3,
     &           ' generations (no termination)')
         end if

C        Selection method
         if (lsel) then
            if (stype .eq. 1) then
               write(*,'(''Using roulette wheel selection'')')
            else
               write(*,
     & '(''Using tournament selection, tournament size = '',i2)')
     &              tsize
            end if
         else
c        Default selection routine
            stype = 1
            write(*,'(''Using DEFAULT of roulette wheel selection'')')
         end if

C        Fitness function
         if (lfit) then
            if (ftype .eq. 1) then
               write(*,15)'linear'
            else if (ftype .eq. 2) then
               write(*,15)'power'
            else if (ftype .eq. 3) then
               write(*,15)'exponential'
            else if (ftype .eq. 4) then
               write(*,15)'tanh'
            else if (ftype .eq. 5) then
               write(*,15)'cos2'
            end if
 15         format('Fitness calculated by ',a,' relationship')  
         else
C        Default fitness function
            ftype = 4
          write(*,'(''Using '',a,'' fitness relationship by DEFAULT'')')
     &           'tanh'
         end if

C        Mating function
         if (lmate) then
            if (mtype .eq. 1) then
               write(*,
     &         '(''Using 1 point crossover about a random point'')')
            else if (mtype .eq. 2) then
               write(*,16)
 16            format('Using 1 point crossover about a point calculated'
     &              ,' from fitness values')
            else if (mtype .eq. 3) then
               write(*,
     &              '(''Using 2 point crossover about random points'')')
            else if (mtype .eq. 4) then
               write(*,17)
 17            format('Using 1 point crossover about random point'
     &              ,' and separating cluster portions')
            else if (mtype .eq. 5) then
               write(*,
     &              '(''Mating cluster using geometric average'')')
            end if
         else
            if (nels .eq. 1) then
c           Default mating routine for single element clusters
               mtype = 3
               write(*,
     &  '(''Using DEFAULT 2 point crossover'')')
            else
c           Default mating routine for dual element clusters
               mtype = 2
               write(*,
     &  '(''Using DEFAULT 1 point weighted crossover'')')
            end if
         end if

C        Mutation method
         if (lmut) then
            if (muttype .eq. 1) then
              write(*,'(''Mutating by cluster replacement'')')
            else if (muttype .eq. 2) then
               write(*,'(''Mutating by half cluster rotation'')')
            else if (muttype .eq. 3) then
               write(*,'(''Mutating by random atomic moves'')')
            else if (muttype .eq. 4) then
               write(*,'(''Mutating by exchange of unlike elements'')')
               if (.not. lmeswaps) then 
                  meswaps = natoms/3
                  write(*,18)meswaps
 18               format('Using DEFAULT mutation exchange of ',i3,
     &                   ' atom swaps') 
               else
                  write(*,19)meswaps
 19               format('Using mutation exchange of ',i3,
     &                   ' atom swaps')
               end if
            end if
         else
c        Default mutation type
            muttype = 1
          write(*,'(''Using DEFAULT mutation by cluster replacement'')')
         end if

C        Check if we are using mutate exchange as well...
         if (lmerate) then
            if (muttype .ne. 4) then
              write(*,'(''Using mutation exchange coupled in...'')')

              if ( .not. lmeswaps) then
C        Set defaults if not assigned yet
                meswaps = natoms/3
                write(*,18)meswaps
              else
                write(*,19)meswaps
              end if
C        Print probability of using mutate exchange within mutations
              write(*,22)merate
 22          format('Using mutate exhange rate probability of ',f4.2,'')
            end if
         else
             merate = 0.d0  
         end if

C	 Include High Energy Mutants
         if (highemut) write(*,
     &        '(''Adding high energy mutants to population'')')
            
C        Include Removal of clusters with similar energies
         if (remove) write(*,
     &     '(''Removing clusters of similar energy from population'')')

c        Stop if cluster energy equals input value
         if (lgmin) then
            write(*,
     & '(''Stopping GA run if lowest energy cluster = '',f13.6)')
     &      gmin
         end if

c        Include predator
         if (pred) then
            write(*,'(a,f13.6)')
     &  'Predator operator will consume clusters with energy = ',epred
         end if

         if (lfile) then
            write(*,'(''Writing output to files named : '',a)')
     &           fname(1:flen)
         else
c        Default filename
            flen = 8
            fname = 'genet000'
       write(*,'(''Writing output to DEFAULT named files : genet000'')') 
         end if
         
         if (write_clust) write(*,
     &        '(''Writing all clusters to file named : '',a)')
     &        fname(1:flen)//'.clusters'

         if (write_energy) write(*,
     &        '(''Writing all clusters to file named : '',a,a)')
     &        fname(1:flen)//'.energies'

         if (write_stats) write(*,'(a,a)')
     &        'Writing out information about new clusters entering',
     &        ' the population'

      end if

      return
      end
      
C     
C     Subroutine to check if variable has already been defined
C
      subroutine test(flag,error,ln)
      implicit none

      logical flag,error
      integer ln

      if (flag) then
         write(*,
     &  '(''ERROR on line '',i2,'', parameter already defined'')')
     &  ln           
         error = .true.
      end if

      return
      end

C
C     Subroutine to read int from tokens
C
      subroutine expect_int(token,ttok,tok,buffer,blen,bloc,lno,
     &     number,error)
      implicit none

C     Character Arrays
      character*(*) token
      character*(*) buffer
C     Integers
      integer ttok
      integer tok
      integer number
      integer lent ! String length
      integer ln
      integer chartoint
C     Integer Arrays
      integer blen(ttok)
      integer bloc(ttok)
      integer lno(ttok)
C     Booleans
      logical error
C     Import function
      external chartoint
      
C     Get next token
      call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok+1)

C    check for equals sign
      if ( token(1:lent) .eq. '=' ) then
C      if ( token .eq. '=' ) then

C        Get next token if it is an equals
         call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok+2)

         if (lent .gt. 0) then
            number = chartoint(token(1:lent),lent)
C             number = chartoint(token,lent)
         else

C     (Original comment)
c     dodgyness - end of file reached - get previous line number

            call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,
     &           tok+1)
            number = -1
         end if

C     (Original comment)            
c     Check for number - value of -1 results from error in chartoint
         if ( number .eq. -1 ) then
            write(*,
     & '(''ERROR on line '',i2,'', number expected after equals'')')ln
            error = .true.
 
C     (Original comment)           
C     return value of tok for next token
            tok = tok+2

         else
C     everything was correct
            tok = tok+3
         end if
      else
         error = .true.
C        (Original comment)
c        check to see if number follows first token
         call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok+1)
      
         if (lent .eq. 0) then
C        (Original comment)
c        EOF reached need to return line number for previous token
            ln = lno(tok)
         end if

         number = chartoint(token(1:lent),lent)
C         number = chartoint(token,lent) 
        
         if ( number .ne. -1) then
c        Missing equals sign
            write(*,
     & '(''ERROR on line '',i2,'', equals sign expected'')')ln
                 
c     return tok
            tok = tok+2
         else
c     no number or equals
            write(*,
     & '(''ERROR on line '',i2,'' no number given after statement '')')
     &           ln       
            tok = tok+1
         end if
      end if
      return
      end

C
C     Subroutine to get double value
C
      subroutine expect_dble(token,ttok,tok,buffer,blen,bloc,lno,dnum,
     &     error)
      implicit none

C     Character Arrays
      character*(*) token
      character*(*) buffer
C     Integers
      integer ttok
      integer tok
      integer number
      integer lent ! String length
      integer ln
C     Integer Arrays
      integer blen(ttok)
      integer bloc(ttok)
      integer lno(ttok)
C     Doubles
      real*8 dnum
      real*8 chartodb
C     Booleans
      logical error
C     Import function
      external chartodb

C     Get next token
      call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok+1)

C    check for equals sign
      if ( token(1:lent) .eq. '=' ) then
C      if ( token .eq. '=' ) then

C        Get next token if it is an equals
         call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok+2)

         if (lent .gt. 0) then
            dnum = chartodb(token(1:lent),lent)
C             dnum = chartodb(token,lent)
         else

C     (Original comment)
c     dodgyness - end of file reached - get previous line number

            call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,
     &           tok+1)
            dnum = -1d0
         end if

C     (Original comment)            
c     Check for number - value of -1 results from error in chartoint
         if ( dnum .eq. -1d0 ) then
            write(*,
     & '(''ERROR on line '',i2,'', number expected after equals'')')ln
            error = .true.

C     (Original comment)           
C     return value of tok for next token
            tok = tok+2

         else
C     everything was correct
            tok = tok+3
         end if
      else
         error = .true.
C        (Original comment)
c        check to see if number follows first token
         call get_token(ttok,buffer,blen,bloc,lno,token,lent,ln,tok+1)

         if (lent .eq. 0) then
C        (Original comment)
c        EOF reached need to return line number for previous token
            ln = lno(tok)
         end if

         dnum  = chartodb(token(1:lent),lent)
C         dnum = chartodb(token,lent)

         if ( dnum .ne. -1d0) then
c        Missing equals sign
            write(*,
     & '(''ERROR on line '',i2,'', equals sign expected'')')ln

c     return tok
            tok = tok+2
         else
c     no number or equals
            write(*,
     & '(''ERROR on line '',i2,'' no number given after statement '')')
     &      ln
            tok = tok+1
         end if
      end if
      return
      end
