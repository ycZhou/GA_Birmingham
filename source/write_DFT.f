c    Genetic Algorithm Program
c    Develoepd by the Birmingham University Cluster Group, Prof.
c    Roy Johnston
c    Commenting added by Andrew Logsdail, Oct 2010

C    This is our subroutine to create the DFT input files. SH 01/11


          SUBROUTINE write_DFT(natoms,na,nb,off,qoff,i,j,
     &                         namea,nameb,
     &                         mass_a,mass_b,nspec,
     &                         prefix,size_cell,
     &                         range,r_ab)
          use hectorstuff
          use DFT
          implicit none

          integer natoms           ! Number of atoms
          integer na               ! Number of atom a
          integer nb               ! Number of atom b
          integer nspec            ! Number of elements
          integer i                ! Count generations or see j
          integer j                ! Count cluster or zero for first Pop
          integer k

          real*8 off(3*natoms+1)   ! Offspring
          real*8 qoff(natoms)      ! Charge (element identifier)
          real*8 bohr(3*natoms)    ! Coordinates in a0 + in cell center
          real*8 mass_a            ! Mass of element a
          real*8 mass_b            ! Mass of element b
          real*8 size_cell         ! Size of primitive unit cell
          real*8 cubic
          real*8 range             ! Used by GPAW DFT call to define cell limits
          real*8 r_ab

          character*60 prefix      ! Prefix of filename
          character*5 na_char      ! Character form of na
          character*5 nb_char      ! Character form of nb
          character*5 i_char       ! Character form of i
          character*5 j_char       ! Character form of j
          character*2 namea        ! Name of element a
          character*2 nameb        ! Name of element b


          logical fexist

c    Converting the integers to string for the filename           

          write(na_char,'(i4)')na
          write(nb_char,'(i4)')nb
          write(i_char,'(i4)')i
          write(j_char,'(i4)')j

C Moved to write_QE area of calculation
c    Convert to mass center coordinate system and place cluster 
c    in center of cubic unit cell with adequate size
C
C          cubic = 0.d0
C
C          call cell_size(off,natoms,mass_a,mass_b,na,nb,qoff,
C     &                   bohr,size_cell)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c    Define filename for monometallic of bimetallic clusters

          if (nb .eq. 0) then
             prefix = trim(namea)//trim(adjustl(na_char))//'_'//
     &                trim(adjustl(i_char))// '_'//
     &                trim(adjustl(j_char))
          else 
             prefix = trim(namea)//trim(adjustl(na_char))//
     &               trim(nameb)//
     &               trim(adjustl(nb_char))//'_'//
     &               trim(adjustl(i_char))//'_'//
     &               trim(adjustl(j_char)) 
          end if    

c    Check if file exists and write to file

          inquire(file = trim(prefix)//'.in',exist=fexist)


          if (dft_type .eq. 1) then

          cubic = 0.d0

          call cell_size(off,natoms,mass_a,mass_b,na,nb,qoff,
     &                   bohr,size_cell)

          call write_QE(natoms,bohr,qoff,namea,nameb,
     &                        mass_a,mass_b,nspec,
     &                        prefix,size_cell,cubic)

          else if (dft_type .eq. 2) then

          !! CREATE GPAW INPUT ON PIPE 23
          call write_GPAW(natoms,off,qoff,
     &                        functional_a,functional_b,namea,nameb,
     &                        scratch,conv_DFT,nspec,
     &                        sys_DFT,ele_DFT,prefix,size_cell,
     &                        openshell,range,r_ab)

          else if (dft_type .eq. 3) then

          !! CREATE G09 INPUT ON PIPE 23
          call write_G09(natoms,off,qoff,namea,nameb,nspec,
     &                        sys_DFT,prefix)


          else if (dft_type .eq. 4) then
 
          !! CREATE NW INPUT ON PIPE 23
          call write_NW(natoms,qoff,off,prefix)


          else if (dft_type .eq. 5) then

          cubic = 0.d0

          call cell_size(off,natoms,mass_a,mass_b,na,nb,qoff,
     &                     bohr,size_cell)

          call write_VASP(natoms,bohr,qoff,namea,nameb,
     &                        mass_a,mass_b,nspec,
     &                        prefix,size_cell,cubic,na_char,
     &                        nb_char)


          end if


          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          SUBROUTINE write_G09(natoms,off,qoff,namea,nameb,nspec,
     &                        sys_DFT,prefix)

          implicit none

          integer natoms           ! Number of atoms
          integer nspec            ! Number of elements
          integer k

          real*8 off(3*natoms+1)   ! Offspring
          real*8 qoff(natoms)      ! Charge (element identifier)

          character*60 prefix      ! Prefix of filename
          character*20 sys_DFT(5)  ! System variables DFT *CURRENTLY UNUSED BUT SUITABLE FOR SOFTCODING OF VARIABLES?* 
                                   ! (1) E_cut_off (2) occupations
                                   ! (3) smearing  (4) degauss
                                   ! (5) number of unpaired electrons
          character*2 namea        ! Name of element a
          character*2 nameb        ! Name of element b

          open(unit=23,status='new',file=trim(prefix)//'.in')

C Inputs
          write(23,'(a)') '%nprocs=8'
          write(23,'(a)') '%mem=13000mb'
C      write(51,'(a)') '%chk=' // trim(name1) // '.chk'
          write(23,'(a)') '# lsda/lanl2mb opt(loose)'
          write(23,'(a)') '# scf(conver=4) int(grid=sg1)'
          write(23,*)
          write(23,'(a)')  trim(prefix)//'.in'
          write(23,*)
          write(23,'(a)') '0 1'

C Atom positions
          do k=1,natoms
            if (qoff(k) .gt. 0.d0) then
              write(23,'(a,1x,3f13.8)') namea,
     &            off(k),
     &            off(k+natoms),
     &            off(k+2*natoms)
            else
              write(23,'(a,1x,3f13.8)') nameb,
     &            off(k),
     &            off(k+natoms),
     &            off(k+2*natoms)
            end if
          end do
          write(23,*)

C Close file
          close(23,status='keep')

          END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          SUBROUTINE write_GPAW(natoms,off,qoff,
     &                        functional_a,functional_b,namea,nameb,
     &                        scratch,conv_DFT,nspec,
     &                        sys_DFT,ele_DFT,prefix,size_cell,
     &                        openshell,range,r_ab)
      use commons, only:nmax,ncmax,maxgen,fnamelen,nspecmax2,maxtoklen,
     &                  maxtokens,tokbufflen
          implicit none

          integer natoms           ! Number of atoms
          integer nspec            ! Number of elements
          integer openshell        ! Is DFT calculation open shell? 
          integer k
          integer l
          integer grid             ! Defines the grid size for GPAW calculation

          real*8 off(3*natoms+1)   ! Offspring
          real*8 qoff(natoms)      ! Charge (element identifier)
C          real*8 bohr(3*natoms)    ! Coordinates in a0 + in cell center
C          real*8 mass_a            ! Mass of element a
C          real*8 mass_b            ! Mass of element b 
          real*8 size_cell         ! Size of primitive unit cell
          real*8 range             ! Breadth of volume of cluster
          real*8 r_ab              ! Mean bond length

          character*60 prefix      ! Prefix of filename
          character*20 functional_a! Functional for element a
          character*20 functional_b! Functional for element b 
          character*20 scratch     ! Specify scratch directory 
          character*20 conv_DFT(3) ! Conv criteria (1) energy (2) force  
                                   ! (3) optimazation step number
          character*20 sys_DFT(5)  ! System variables DFT 
                                   ! (1) E_cut_off (2) occupations
                                   ! (3) smearing  (4) degauss
                                   ! (5) number of unpaired electrons
          character*20 ele_DFT(2)  ! Electron variables DFT
                                   ! (1) SCF convergence (2) mixing SCF
          character*2 namea        ! Name of element a
          character*2 nameb        ! Name of element b
CCCCCCCCCC Characters for GPAW input file CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          character*(nmax*2) atms 
          character*(nmax*2) temp_atms        ! Atoms in input file
          character*20 c_atms(3)              ! Current atom positions
          character*(1+65*nmax) positions     ! Positions in input file
          character*(1+65*nmax) temp_positions! Temp positions for catenating words
          character*65 c_positions            ! Positions of current atom

          atms = ""
          positions = ""
          size_cell =range*r_ab+10.0
          grid =nint((size_cell/0.3)/8.0) * 8

          do k = 1, natoms
C Firstly prepare labelling string
             temp_atms = atms
             if (qoff(k) .eq. 1.d0) then
               atms = trim(temp_atms)//trim(namea)
             else
               atms = trim(temp_atms)//trim(nameb)
             end if

C Now we prepare the positional string at the end
C Firstly lets clear the values
C And then in the same loop give them their values.
 
             do l = 1, 3
                c_atms(l) = ""
                write(c_atms(l),'(f13.8)') off(k+((l-1)*(natoms))) 
             end do            

             c_positions='('//trim(c_atms(1))//','
     &                       //trim(c_atms(2))//','
     &                       //trim(c_atms(3))//')'

             temp_positions = positions
             if (k .eq. natoms) then
                positions=trim(temp_positions)//
     &                    trim(c_positions)
             else
                positions=trim(temp_positions)//
     &                    trim(c_positions)//', '
             end if 
          end do

C Open file for writing
          open(unit=23,status='new',file=trim(prefix)//'.in')

C Prepare input environment
          write(23,'(a)') "#!/usr/bin/env python"
          write(23,'(a)') 'import sys'
          write(23,'(a)') 'from ase import *'
          write(23,'(a)') 'from gpaw import *'
          write(23,'(a)') 'from ase.optimize.bfgslinesearch '//
     &                      'import BFGSLineSearch'
          write(23,*)
          write(23,'(a,f10.6),') 'a =',size_cell
          write(23,*)

CCCCCCCCC INSERT METHOD TO WRITE ALL ATOMS HERE

          write(23,'(a)') "mol = Atoms('"//trim(atms)//"',"
          write(23,'(a)') "            positions=["//
     &                          trim(positions)//"])"

CCCCCCCCC INSERT METHOD TO WRITE ALL ATOMS HERE

          write(23,'(a)') 'mol.set_cell((a,a,a))'
          write(23,'(a)') 'mol.set_pbc((False,False,False))'
          write(23,'(a)') 'mol.center()'
          write(23,*)
          !! We have prepared the molecule 

          !! Prepare calculator
          write(23,'(a)') "calc = GPAW(mode='lcao', basis='"
     &         //trim(scratch)//"',"//
     &         "poissonsolver=PoissonSolver(relax='GS',eps=1e-7),"
          write(23,'(a,i3,a,i3,a,i3,a)') '        gpts=(',
     &         grid,',',grid,',',grid,'),'
          write(23,'(a)') "        xc='"//trim(functional_a)//"',"

CCCCCCCCCCCCCC
          !! Paired system
          if (openshell .eq. 0) then
             write(23,'(a)') '        spinpol=False,'
             write(23,'(a)') '        mixer=Mixer(beta='//
     &                                trim(conv_DFT(1))//
     &                                ', nmaxold='//
     &                                trim(conv_DFT(2))//
     &                                ', weight='//
     &                                trim(conv_DFT(3))//'),'

          else
          !! Unpaired system
             write(23,'(a)') '        spinpol=True,'
             write(23,'(a)') '        mixer=MixerSum(beta='//
     &                                trim(conv_DFT(1))//
     &                                ', nmaxold='//
     &                                trim(conv_DFT(2))//
     &                                ', weight='//
     &                                trim(conv_DFT(3))//'),'
          end if
CCCCCCCCCCCCCC

          write(23,'(a)') '        maxiter = '//trim(sys_DFT(1))//','
          write(23,'(a)') '        occupations=FermiDirac(width='//
     &                             trim(sys_DFT(2))//'),'
          write(23,'(a)') '        nbands='//trim(sys_DFT(3))//','
          write(23,'(a)') "        eigensolver='"//trim(sys_DFT(4))
     &                             //"',"
          write(23,'(a)') "        convergence = {'energy': "//
     &                             trim(ele_DFT(1))//", 'density': "//
     &                             trim(ele_DFT(2))//"},"
          write(23,'(a)') "        txt='"//trim(prefix)//".out')"

          !! Now all we have to do is run the minimisation
          write(23,*)
          write(23,'(a)') 'mol.set_calculator(calc)'
          write(23,'(a)') 'opt = BFGSLineSearch(mol)' 
          write(23,'(a)') 'try:'
          write(23,'(a)') '    opt.run(fmax='//
     &                         trim(functional_b)//')'
          write(23,'(a)') 'except:'
          write(23,'(a)') '    sys.exit()'
C          write(23,'(a)') "#write(name+'_final.xyz', mol)"
    
          close(23,status='keep')

          END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          SUBROUTINE write_QE(natoms,bohr,qoff,namea,nameb,
     &                        mass_a,mass_b,nspec,
     &                        prefix,size_cell,cubic)
          use HECTORSTUFF
          use commons, only: nels
          implicit none
 
          integer natoms           ! Number of atoms
          integer nspec            ! Number of elements
          integer k

          real*8 qoff(natoms)      ! Charge (element identifier)
          real*8 bohr(3*natoms)    ! Coordinates in a0 + in cell center
          real*8 mass_a            ! Mass of element a
          real*8 mass_b            ! Mass of element b
          real*8 size_cell         ! Size of primitive unit cell
          real*8 cubic

          character*60 prefix      ! Prefix of filename
          character*2 namea        ! Name of element a
          character*2 nameb        ! Name of element b
          character*50 line

          open(unit=24,file='DFT.in')
          open(unit=23,file=trim(prefix)//'.in')

c Read in the DFT.in file and print to new input file.

2000      read(24,'(A50)',end=1000) line

c Print out the unit cell for each input generated.

          if(trim(line)=="BCGAcell") then
            write(23,'(f5.2,1x,f5.2,1x,f5.2)')size_cell,cubic,cubic
            write(23,'(f5.2,1x,f5.2,1x,f5.2)')cubic,size_cell,cubic
            write(23,'(f5.2,1x,f5.2,1x,f5.2)')cubic,cubic,size_cell

c Print out coordinates for each input generated.

          else if(trim(line)=="BCGAcoords") then
            do k = 1, natoms
            if (qoff(k) .eq. 1.d0 .or. nels .eq. 1) then
            write(23,'(a2,1x,3f13.8)')trim(adjustl(namea)),bohr(k),
     &                  bohr(k+natoms),bohr(k+2*natoms)
            else
            write(23,'(a2,1x,3f13.8)')trim(adjustl(nameb)),bohr(k),
     &                  bohr(k+natoms),bohr(k+2*natoms)
            end if
          end do

c Finish printing DFT.in to new input file.

          else
            write(23,'(a)') line
          endif
          goto 2000
1000      close(23)
          close(24)


          END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE write_NW(natoms,qoff,off, prefix)
      use potentials, only: nspec, namea,nameb
      use commons, only: nels
      implicit none

      integer natoms           ! Number of atoms
      integer k
      real*8 off(3*natoms+1)   ! Offspring
      real*8 qoff(natoms)      ! Charge (element identifier)
      
      character*60 prefix      ! Prefix of filename
      character*50 line

      open(unit=24,file='DFT.nw')
      open(unit=23,file=trim(prefix)//'.nw')

 2000 read(24,'(A50)',end=1000) line
      if(trim(line)=="BCGAcoords") then
         do k= 1, natoms
            if (qoff(k) .eq. 1.d0 .or. nels.eq.1) then
               write(23,'(a2,1x,3f13.8)')trim(adjustl(namea)),off(k),
     &                off(k+natoms),off(k+2*natoms)
            else
               write(23,'(a2,1x,3f13.8)')trim(adjustl(nameb)),off(k),
     &               off(k+natoms),off(k+2*natoms)
            end if
         end do
      else
         write(23,'(a)') line
      endif
      goto 2000
 1000 close(23)
      close(24)
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          SUBROUTINE write_VASP(natoms,bohr,qoff,namea,nameb,
     &                        mass_a,mass_b,nspec,
     &                        prefix,size_cell,cubic,na_char,
     &                        nb_char)
          use HECTORSTUFF
          use commons, only: nels
          implicit none

          integer natoms           ! Number of atoms
          integer nspec            ! Number of elements
          integer k

          real*8 qoff(natoms)      ! Charge (element identifier)
          real*8 bohr(3*natoms)    ! Coordinates in a0 + in cell center
          real*8 mass_a            ! Mass of element a
          real*8 mass_b            ! Mass of element b
          real*8 size_cell         ! Size of primitive unit cell
          real*8 cubic

          character*60 prefix      ! Prefix of filename
          character*2 namea        ! Name of element a
          character*2 nameb        ! Name of element b
          character*5 na_char
          character*5 nb_char
          character*50 line

          open(unit=24,file='DFT.in')
          open(unit=23,file='POSCAR')

          write(23,'(a)')trim(prefix)

c Read in the POSCAR file and print to new input file.

2000      read(24,'(A50)',end=1000) line

c Print out the unit cell for each input generated.

          if(trim(line)=="BCGAcell") then
            write(23,'(f5.2,1x,f5.2,1x,f5.2)')size_cell,cubic,cubic
            write(23,'(f5.2,1x,f5.2,1x,f5.2)')cubic,size_cell,cubic
            write(23,'(f5.2,1x,f5.2,1x,f5.2)')cubic,cubic,size_cell
            if (trim(nb_char) .eq. "0") then
               write(23,'(a)')trim(na_char)
            else
               write(23,'(a)')trim(na_char)//' '//trim(nb_char)
            endif

c Print out coordinates for each input generated.

          else if(trim(line)=="BCGAcoords") then
            do k = 1, natoms
            if (qoff(k) .eq. 1.d0 .or. nels .eq. 1) then
            write(23,'(1x,3f13.8)')bohr(k),
     &                  bohr(k+natoms),bohr(k+2*natoms)
            endif
            if (qoff(k) .ne. 1.d0) then
            write(23,'(1x,3f13.8)')bohr(k),
     &                  bohr(k+natoms),bohr(k+2*natoms)
            end if
          end do

c Finish printing POSCAR to new input file.

          else
            write(23,'(a)') line
          endif
          goto 2000
1000      close(23)
          close(24)

          CALL system( 'mkdir '//trim(prefix) )
          CALL system( 'cp INCAR '//trim(prefix)//'/' )
          CALL system( 'mv POSCAR '//trim(prefix)//'/' )
          CALL system( 'cp POTCAR '//trim(prefix)//'/' )
          CALL system( 'cp KPOINTS '//trim(prefix)//'/' )

          END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          SUBROUTINE cell_size(off,natoms,mass_a,mass_b,na,nb,qoff,
     &                          bohr,size_cell)

          use DFT
          use commons, only: nels
          implicit none

          integer natoms
          integer na
          integer nb
          integer i

          real*8 off(3*natoms+1)
          real*8 qoff(natoms)
          real*8 bohr(3*natoms)
          real*8 mass_a
          real*8 mass_b
          real*8 Mass              ! Total mass of cluster
          real*8 max_bohr
          real*8 a0
          real*8 cell_para
          real*8 size_cell
          real*8 R(3)              ! Center of mass vector

c Set parameters for calculation

          a0 = 0.5291772           ! Bohr to Angstrom
          cell_para = 16.d0/a0     ! distance between clusters in a0
          if (nels.eq.2) then
             Mass = na * mass_a + nb * mass_b
          else
             mass = natoms*mass_a
          endif

          R = 0.d0

c Find center of mass
         
          do i = 1, natoms
             if (qoff(i) .eq. 1.d0 .or. nels .eq. 1) then
               R(1) = R(1)+off(i)*mass_a/Mass
               R(2) = R(2)+off(i+natoms)*mass_a/Mass
               R(3) = R(3)+off(i+2*natoms)*mass_a/Mass
             else
               R(1) = R(1)+off(i)*mass_b/Mass
               R(2) = R(2)+off(i+natoms)*mass_b/Mass
               R(3) = R(3)+off(i+2*natoms)*mass_b/Mass
             end if
          end do


c Convert to Bohr and to center of mass coordinate system
c VASP coordinates in angstrom

          if (dft_type .eq. 5) then
          do i = 1, natoms
             bohr(i) = (off(i)-R(1))
             bohr(i+natoms) = (off(i+natoms)-R(2))
             bohr(i+2*natoms) = (off(i+2*natoms)-R(3))
          end do
          else
          do i = 1, natoms
             bohr(i) = (off(i)-R(1)) / a0
             bohr(i+natoms) = (off(i+natoms)-R(2)) / a0
             bohr(i+2*natoms) = (off(i+2*natoms)-R(3)) / a0
          end do
          endif

c Find cell size

          max_bohr = maxval(abs(bohr))

c VASP cell in angstrom else bohr

          if (dft_type .eq. 5) then
            size_cell = max_bohr + (cell_para * a0)
          else
            size_cell = max_bohr + cell_para
          endif

c Place cluster in the middle of a cubic unit cell

          do i = 1, natoms
             bohr(i) = bohr(i) + size_cell/2.d0
             bohr(i+natoms) = bohr(i+natoms) + size_cell/2.d0
             bohr(i+2*natoms) = bohr(i+2*natoms) + size_cell/2.d0
          end do

c Convert VASP coordinates to direct

          if (dft_type .eq. 5) then
            do i = 1, natoms
            bohr(i) = bohr(i) / (cell_para * a0)
            bohr(i+natoms) = bohr(i+natoms) / (cell_para * a0)
            bohr(i+2*natoms) = bohr(i+2*natoms) / (cell_para * a0)
            end do
          endif

          END
