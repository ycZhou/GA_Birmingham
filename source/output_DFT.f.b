c    Genetic Algorithm Program
c    Develoepd by the Birmingham University Cluster Group, Prof.
c    Roy Johnston

c    This is our subroutine to read the DFT output files and
c    transform it into the format required for the GA. SH 01/11

C    Updates 07/11: AJL
C    Add functionality for GPAW and G09
C    Split everything down into subroutines as I like my C objects

          SUBROUTINE output_DFT(natoms,dft_type,prefix,qoff,off,
     &                          namea,nameb,size_cell,iflag,
     &                          jobDone,atn_a,atn_b)

          implicit none

          integer natoms             ! Number of atoms
          integer dft_type	     ! DFT Type where (1) QE (2) GPAW (3) G09
          integer atn_a              ! Atomic Number of A
          integer atn_b              ! Atomic Number of B

          real*8 off(3*natoms+1)     ! Offspring
          real*8 qoff(natoms)        ! Charge (element identifier)
          real*8 jobDone             ! Is job done ?
          real*8 size_cell
                
          character*60 prefix        ! Prefix of filename
          character*2 namea          ! Name of element a
          character*2 nameb          ! Name of element b
          character*100 cwd          !current working directory at point of reading output file


          logical fexist 
          logical iflag              ! Indicates the converged DFT-Calc 
          logical file_exists

          jobDone = 0.d0

c          CALL GETCWD(cwd)
c          write (*,*) TRIM(cwd)
c          write (*,*) prefix
c          write (*,*) trim(prefix)//'.out'

          INQUIRE(FILE=trim(prefix)//'.out',EXIST=file_exists) 
          if(file_exists .eqv. .true.) then
c           if (fexist) then
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

            end if

            close(24,status='keep')
          else
            write(*,*) 'Error reading DFT output'
            write(*,*) 'File not found'

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
              read(24,'(a)',end=333)jobstat
              lineno1=index(jobstat,'   JOB DONE.')
              if (lineno1.ne.0)then
                 jobDone = 1.d0
                 loop = .false.
                 write(*,'(a)')'JOBDONE_WAS_FOUND'
              else
                 write(*,'(a)')'JOBDONE WASNT WRITTEN, AS IF'
              end if
          end do

c    Check if Job is Done. If not get new random structure.
c333       if (jobDone .eq. 1.d0) then
333             rewind (24)

c    Read final DFT energy and coordinates. What happens if coordinates
c    are of magnitude of 1e+02 a0? Not sure yet! Should be ok for small
c    clusters.

C            if (jobDone .eq. 1.d0) then

              loop = .true.
              do while (loop)
                  read(24,'(a)',end=555)jobstat
                  lineno1=index(jobstat,
     &                  '     End of BFGS Geometry Optimization')
cHERE I WANT THE COORDS TO BE WRITTEN OUT WHATEVER. IF NOT CONVERGED,
C IT WILL STILL RESTART, BUT I WANT TO SEE THE COORDS, SO I'VE REMOVED
C THE IF LOOP FOR READING COORDS IN
c                  if (lineno1 .ne. 0) then
                   read(24,'(a)')
                   read(24,'(24x,f17.10)')energy_Ry
                   read(24,'(a)')
                   read(24,'(a)')
                   read(24,'(a)')
                   do i =1,natoms
                   read(24,'(A2,4x,f12.9,2x,f12.9,2x,f12.9)')ename(i),
     &             bohr(i),bohr(i+natoms),bohr(i+2*natoms)
                   write(*,'(f12.9,2x,f12.9,2x,f12.9)')bohr(i),
     &             bohr(i+natoms),bohr(i+2*natoms)
               end do
c ONCE READ IN THE COORDS, WRITE TO * THE COORDS 
C BASICALLY THIS MEANS IF JOBDONE WAS THERE, I GET THE COORDS OF THE
C LAST GEOMETRY. ALSO MEANS IF I GET THE GEOMETRY, JOBDONE WORKED, SO IT 
C SHOULD NOT HAVE CANCELLED, AND THE FAULT IS IN THE SECOND GREP BIT
                      do i=1,natoms
                           write(*,'(f12.9,2x,f12.9,2x,f12.9)')bohr(i),
     &                     bohr(i+natoms),bohr(i+2*natoms)
                      enddo
                      if (lineno1 .ne. 0) then
                        indicator = 1
                        iflag = .false.
                        loop = .false.
                      end if
              end do

c    If SCF does not converge find structure before last step
c    and get energy and coordinates. I.E. Not Complete!

555           if (indicator .eq. 0) then

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

c          else
        
C    Return with unknown error to create random new structure
C    jobDone should still be set to zero so we won't reassign
c              iflag = .true.
c          end if

          END

