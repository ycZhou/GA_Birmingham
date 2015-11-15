c    Genetic Algorithm Program
c    Develoepd by the Birmingham University Cluster Group, Prof.
c    Roy Johnston

c    This is our subroutine to read the DFT output files and
c    transform it into the format required for the GA. SH 01/11

C    Updates 07/11: AJL
C    Add functionality for GPAW and G09
C    Split everything down into subroutines as I like my C objects

          SUBROUTINE output_DFT(natoms,prefix,qoff,off,
     &                          namea,nameb,size_cell,iflag,
     &                          jobDone)
          use DFT
          implicit none

          integer natoms             ! Number of atoms

          real*8 off(3*natoms+1)     ! Offspring
          real*8 qoff(natoms)        ! Charge (element identifier)
          real*8 jobDone             ! Is job done ?
          real*8 size_cell
                
          character*60 prefix        ! Prefix of filename
          character*2 namea          ! Name of element a
          character*2 nameb          ! Name of element b

          logical fexist 
          logical iflag              ! Indicates the converged DFT-Calc 
          logical file_exists

          jobDone = 0.d0


c Set up read VASP output

          if (dft_type .eq. 5) then
          INQUIRE(FILE=prefix,EXIST=file_exists)
          if(file_exists .eqv. .true.) then

          open(unit=25,file=trim(prefix)//'/OUTCAR',status='old')
          open(unit=26,file=trim(prefix)//'/CONTCAR',status='old')

          call output_VASP(natoms,prefix,qoff,off,
     &                        namea,nameb,iflag,
     &                        jobDone)

          close(25,status='keep')
          close(26,status='keep')
          else
          write(*,*) 'Error reading DFT output'
          write(*,*) 'File not found'
          endif

c Set up read in other output

          else
c          if (fexist) then
          INQUIRE(FILE=trim(prefix)//'.out',EXIST=file_exists)
          if(file_exists .eqv. .true.) then
                            
            open(unit=24,file=trim(prefix)//'.out',status='old')
            rewind(24)

            if (dft_type .eq. 1) then
  
               call output_QE(natoms,prefix,qoff,off,
     &                        namea,nameb,size_cell,iflag,
     &                        jobDone)

            else if (dft_type .eq. 2) then

               call output_GPAW(natoms,prefix,qoff,off,
     &                          namea,nameb,iflag,jobDone)

            else if (dft_type .eq. 3) then

            !! G09 OUTPUT MANIPULATION
               call output_G09(natoms,prefix,qoff,off,
     &                         namea,nameb,iflag,jobDone,atn_a,atn_b)


            else if (dft_type .eq. 4) then

            !! NW OUTPUT MANIPULATION
               call output_NW(natoms,prefix,qoff,off,
     &                        namea,nameb,iflag,
     &                        jobDone)


            end if

            close(24,status='keep')
          else
            write(*,*) 'Error reading DFT output'
            write(*,*) 'File not found'


          end if
          end if 

          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          SUBROUTINE output_G09(natoms,prefix,qoff,off,
     &                           namea,nameb,iflag,jobDone,atn_a,atn_b)

          integer natoms             ! Number of atoms
          integer i
          integer stat
          integer indicator          ! Indicates the converged DFT-Calc
          integer atn(natoms)        ! atomic numbers from output file
          integer atn_a              ! Atomic Number A
          integer atn_b              ! Atomic Number B

          real*8 off(3*natoms+1)     ! Offspring
          real*8 qoff(natoms)        ! Charge (element identifier)
          real*8 jobDone             ! Is job done ?

          character*70 line          ! Line from output file
          character*60 prefix        ! Prefix of filename
          character*2 namea          ! Name of element a
          character*2 nameb          ! Name of element b
          character*20 temp          ! For filename manipulation

          logical fexist
          logical iflag              ! Indicates the converged DFT-Calc

          !! G09 OUTPUT MANIPULATION
          !! Lifted from work by Shunjiang Yang, ANL, 2011

          ! Set everything to zero
          stat=0
          indicator=0

          do
            read(24,fmt='(A60)',iostat=stat) line
            ! Check if search has finished without convergence
            if (stat /= 0) exit
            if (line(2:41) ==
     &         'Predicted change in Energy=          NaN') exit

            ! Get energy of cluster for completed search
            if (line(2:10) == 'SCF Done:') then
                  temp = line(24:44)
                  read(temp,'(F20.6)') off(3*natoms+1)

            ! Grab coordinates
            else if (line(27:44) == 'Input orientation:') then
            
                  off = 0.d0

                  read(24,*)
                  read(24,*)
                  read(24,*)
                  read(24,*)

                  do i=1,natoms
                    read(24,*) temp,atn(i),temp,
     &                  off(i),
     &                  off(i+natoms),
     &                  off(i+2*natoms)
                  end do

                  jobDone = 1.d0
                  indicator = 1
                  iflag = .false.

c    Write charges again. Not sure if DFT code rearranges atoms while
c    running. Perhaps not necessary.
                  qoff = 0.d0

                  do i=1,natoms
                      if (atn(i).eq.atn_a)then
                          qoff(i)=1.d0
                      else
                          qoff(i)=2.d0
                      end if
                  end do
            end if
          end do

          ! Raise error if not coordinates found
          if (indicator .eq. 0) iflag = .true.

          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          SUBROUTINE output_GPAW(natoms,prefix,qoff,off,
     &                           namea,nameb,iflag,jobDone)

          integer natoms             ! Number of atoms
          integer i
          integer stat
          integer indicator          ! Indicates the converged DFT-Calc

          real*8 off(3*natoms+1)     ! Offspring
          real*8 qoff(natoms)        ! Charge (element identifier)
          real*8 jobDone             ! Is job done ?

          character*70 line	     ! Line from output file
          character*60 prefix        ! Prefix of filename
          character*2 ename(natoms)  ! Name of element
          character*2 namea          ! Name of element a
          character*2 nameb          ! Name of element b
          character*20 temp	     ! For filename manipulation

          logical fexist
          logical iflag              ! Indicates the converged DFT-Calc

          !! GPAW OUTPUT MANIPULATION
          !! Lifted from work by Shunjiang Yang, ANL, 2011

          ! Assign coordinates and q
          stat = 0
          indicator = 0

          do
            read(24,fmt='(A70)',iostat=stat) line
            ! End of search
            if (stat /= 0) exit

            ! Energy of cluster
            if (line(1:12) == 'Free Energy:') then
              temp = line(13:25)
              read(temp,'(F20.6)') off(3*natoms+1)

            ! Grab coordinates
            else if (line(1:10) == 'Positions:') then

              off = 0.d0

              do i=1,natoms
                read(24,*) temp,ename(i),off(i),
     &          off(i+natoms),off(i+2*natoms)
              end do

              jobDone = 1.d0
              indicator = 1
              iflag = .false.

c    Write charges again. Not sure if DFT code rearranges atoms while
c    running. Perhaps not necessary.
              qoff = 0.d0
          
              do i=1,natoms
                  if (ename(i).eq.namea)then
                      qoff(i)=1.d0
                  else
                      qoff(i)=2.d0
                  end if
              end do
            end if
          end do

          ! Raise error if not coordinates found
          if (indicator .eq. 0) iflag = .true.

          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          SUBROUTINE output_QE(natoms,prefix,qoff,off,
     &                          namea,nameb,size_cell,iflag,
     &                          jobDone)

          implicit none

          integer natoms             ! Number of atoms
          integer dft_type           ! DFT Type where (1) QE (2) GPAW (3) G09
          integer lineno1            ! Integer of keyword indicator
          integer i
          integer indicator          ! Indicates the converged DFT-Calc

          real*8 off(3*natoms+1)     ! Offspring
          real*8 qoff(natoms)        ! Charge (element identifier)
          real*8 bohr(3*natoms)      ! Coordinates in a0 + in cell center
          real*8 energy_Ry           ! Temp energy in Ry
          real*8 a0
          real*8 Ry_eV
          real*8 size_cell
          real*8 jobDone             ! Is job done ?
                
          character*60 prefix        ! Prefix of filename
          character*50 jobstat       ! Keyword indicator
          character*2 ename(natoms)  ! Name of element
          character*2 namea          ! Name of element a
          character*2 nameb          ! Name of element b

          logical fexist 
          logical iflag              ! Indicates the converged DFT-Calc 
          logical loop               ! We are going to use this to remove the goto's

          a0 = 0.5291772             ! Bohr to Angstrom
          Ry_eV = 13.6056923         ! Rydberg to eV 
          bohr = 0.d0
          indicator = 0

c    Read output file until keyword is found.
          loop = .true.
          do while (loop)
              read(24,'(a)',end=222)jobstat
              lineno1=index(jobstat,'   JOB DONE.')
              if (lineno1.ne.0)then
                 jobDone = 1.d0
                 loop = .false.
              write(*,*) 'job is done'
              end if
          end do


222       rewind(24)
          loop = .true.
          do while (loop)
              read(24,'(a)',end=333)jobstat
              lineno1=index(jobstat,'     The maximum number')
              if (lineno1.ne.0)then
                 jobDone = 0.d0
                 loop = .false.
              write(*,*) 'job is not done'
              end if
          end do



c    Check if Job is Done. If not get new random structure.
333       if (jobDone .eq. 1.d0) then
              rewind (24)
               write(*,*) 'yup, still done'
c    Read final DFT energy and coordinates. What happens if coordinates
c    are of magnitude of 1e+02 a0? Not sure yet! Should be ok for small
c    clusters.

C            if (jobDone .eq. 1.d0) then

              loop = .true.
              do while (loop)
                  read(24,'(a)',end=555)jobstat
                  lineno1=index(jobstat,
     &                  '     End of BFGS Geometry Optimization')
                  if (lineno1 .ne. 0) then
                      write(*,*) 'end of bfgs is there. good.'
                      read(24,'(a)')
                      read(24,'(24x,f17.10)')energy_Ry
                      read(24,'(a)')
                      read(24,'(a)')
                      read(24,'(a)')
                      do i =1,natoms
                      read(24,'(a,4x,f12.9,2x,f12.9,2x,f12.9)')ename(i),
     &                     bohr(i),bohr(i+natoms),bohr(i+2*natoms)
                      end do
                      indicator = 1
                      iflag = .false.
                      loop = .false.
                      write(*,*) ' indicator should be 1 and is' 
                      write(*,*) indicator
                      write(*,*) 'iflag should be false and is' 
                      write(*,*) iflag
                  end if
              end do

c    If SCF does not converge find structure before last step
c    and get energy and coordinates. I.E. Not Complete!

555           if (indicator .eq. 0) then
c              write(*,*) 'indicator is 0 loop says'
c              write(*,*) indicator
              
C Identify this as an error so we can use this geometry at restart
                  iflag = .true.
              end if

c    Translating cluster back to center of mass coordinate system.
              bohr=bohr-size_cell/2.d0

c    Convert to Angstrom/eV and write in off
              off = 0.d0

              do i=1,natoms
                  off(i)=bohr(i)*a0
                  off(i+natoms)=bohr(i+natoms)*a0
                  off(i+2*natoms)=bohr(i+2*natoms)*a0
              end do

              off(3*natoms+1)=energy_Ry*Ry_eV

c    Write charges again. Not sure if DFT code rearranges atoms while
c    running. Perhaps not necessary.
              qoff = 0.d0

              do i=1,natoms
                  if (ename(i).eq.namea)then
                      qoff(i)=1.d0
                  else
                      qoff(i)=2.d0
                  end if
              end do

          else
        
C    Return with unknown error to create random new structure
C    jobDone should still be set to zero so we won't reassign
              iflag = .true.
              write(*,*) 'not converged, leaving output_DFT.f'
       end if
           write(*,*)'at end of outputDFT, indicator is'
           write(*,*) indicator
           write(*,*) 'at end of outputDFT, iflag is'
           write(*,*) iflag
          END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         SUBROUTINE output_NW(natoms,prefix,qoff,off,
     &                        namea,nameb,iflag,
     &                        jobDone)
        
         implicit none

         integer natoms             ! Number of atoms
         integer lineno1            ! Integer of keyword indicator
         integer i
         integer indicator          ! Indicates the converged DFT-Calc
         integer stat

         real*8 off(3*natoms+1)     ! Offspring
         real*8 qoff(natoms)        ! Charge (element identifier)
         real*8 energy              ! Energy in eV
         real*8 jobDone             ! Is job done ?

         logical iflag              ! Indicates the converged DFT-Calc 
         logical loop               ! We are going to use this to remove the goto's
         character*60 prefix        ! Prefix of filename
         character*50 jobstat       ! Keyword indicator
         character*2 ename(natoms)  ! Name of elements
         character*2 namea          ! Name of element a
         character*2 nameb          ! Name of element b

         indicator = 0
         stat = 0

!    Read output file until keyword is found
         loop = .true.
         do while (loop)
             read(24,'(a)',iostat=stat)jobstat
             if (stat .ne. 0) exit
             lineno1=index(jobstat,'      Optimization converged')
             if (lineno1.ne.0)then
                jobDone = 1.d0
                do i= 1,5
                  read(24,'(a)')
                end do
                read(24,'(7x,f17.10)')off(3*natoms+1)
                do i= 1,11
                   read(24,'(a)')
                end do
                do i=1,natoms
                read(24,'(6x,a,26x,f14.9,x,f14.9,x,f14.9)')ename(i),
     &               off(i),off(i+natoms),off(i+2*natoms)
                end do
                iflag = .false.
                indicator = 1
                loop = .false.
             write(*,*) 'job is done'
             end if
         end do

         if (indicator .eq. 0) then
              iflag=.true.
              write(*,*)'job failed. resubmit new job'
         end if

         END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         SUBROUTINE output_VASP(natoms,prefix,qoff,off,
     &                        namea,nameb,iflag,
     &                        jobDone)

         implicit none

         integer natoms             ! Number of atoms
         integer dft_type           ! DFT Type where (1) QE (2) GPAW (3) G09
         integer lineno1            ! Integer of keyword indicator
         integer i
         integer indicator          ! Indicates the converged DFT-Calc
         integer stat

         real*8 off(3*natoms+1)     ! Offspring
         real*8 qoff(natoms)        ! Charge (element identifier)
         real*8 bohr(3*natoms)      ! Coordinates in a0 + in cell center
         real*8 energy              ! Temp energy in Ry
         real*8 a0
         real*8 Ry_eV
         real*8 size_cell
         real*8 jobDone             ! Is job done ?

         character*60 prefix        ! Prefix of filename
         character*100 jobstat       ! Keyword indicator
         character*2 ename(natoms)  ! Name of element
         character*2 namea          ! Name of element a
         character*2 nameb          ! Name of element b

         logical fexist
         logical iflag              ! Indicates the converged DFT-Calc
         logical loop

         indicator = 0
         stat = 0

         write(*,*)'INSIDE OUTPUT DFT HERE'
         
         ! Reads OUTCAR until it string is found
         loop = .TRUE. 
         DO WHILE (loop)
                READ(25,'(A)', END = 10) jobstat
                lineno1=INDEX(jobstat,'reached required accuracy -
     &           stopping')
                        IF (lineno1 .NE. 0) THEN
                                jobDone = 1.d0
                                loop = .false. 
                                WRITE(*,*) 'JOB IS DONE'
                        END IF 
         END DO 

10       REWIND(25)

c Need to add search for non-converged calc 

c Searching for final energy
         loop = .TRUE.
         DO WHILE (loop)
                READ(25,'(A)', END = 10) jobstat 
                lineno1=INDEX(jobstat,'  energy  without entropy=')
                        IF (lineno1 .NE. 0) THEN 
                                WRITE(*,*) 'Looking for the energy' 
                                READ(jobstat(65:83),'(F16.6)') energy
                                WRITE(*,*) 'Found it' 
c                                loop = .FALSE.
                        ELSE 
c                                WRITE(*,*) 'Energy not found'
                        END IF
         END DO 

         loop = .TRUE.
         DO WHILE (loop)
                READ(26,'(A)')
                READ(26,'(A)')
                READ(26,'(A)')
                READ(26,'(A)')
                READ(26,'(A)')
                READ(26,'(A)')
                READ(26,'(A)')
                READ(26,'(A)')
                DO i = 1, natoms
                        READ(26,'(2x,F18.16,2x,F18.16,2x,F18.16)')
     &                  off(i), off(i*natoms), off(i+2*natoms)
                END DO 
                iflag = .false.
                indicator = 1 
                loop = .false. 
                WRITE(*,*) 'Job is done'
         END DO 

         IF ( indicator .EQ. 0 ) THEN 
                iflag = .true. 
                WRITE(*,*) 'JOB HAS FAILED, RESUBMIT NEW JOB'
         END IF 

c Read output for final xyz coords
c           if (jobDone .eq. 1.d0) then
c           loop = .true.
c           do while(loop)
c Open CONTCAR file
C             read(26,'(a)',iostat=stat)jobstat
C             do i=1,8
C               read(26,'(a)')
C             end do
C             do i =1,natoms
C             read(26,'(2x,f18.16,2x,f18.16,2x,f18.16)')bohr(i),
C     &                     bohr(i+natoms),bohr(i+2*natoms)
C             end do
C             iflag = .false.
C             indicator = 1
C             loop = .false.
C             write(*,*) 'job is done'
C             end do
C           end if
                
c PJ Search for Energy
C           if (jobDone .eq. 1.d0) then
C            rewind (25)
C             loop = .true.
C             do while(loop)
C               read(25,'(a)',iostat=stat)jobstat
C               if (stat .ne. 0) exit
C               lineno1=index(jobstat,'  energy  without entropy=')
C               if (lineno1.ne.0)then
c May need adjusting
C                 read(jobstat(63:78),'(F16.6)') energy
C               end if
C             end do
C            end if

C          write(*,*)'ENERGY NOT FOUND'
         END
