c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Subrotuine to read in DFT parameters SH 01/11

      MODULE hectorstuff
      implicit none 
      integer :: nproc=0
      integer :: ntask=0
      character*20 :: machine
      END MODULE


      SUBROUTINE read_DFT(inputfile)
      USE hectorstuff
      use DFT
      implicit none
      character*(*) inputfile   ! Input filename


      ! Initialise so we don't get any odd errors
      ! Initialise complete

      open(unit=22,file=trim(adjustl(inputfile)),status='old')

c     Read number of elements
      
      read(22,*)dft_program
     
      if (dft_program .eq. 'QE') then

        dft_type = 1
        call read_QE()

      else if (dft_program .eq. 'GPAW') then
      
        dft_type = 2
        call read_GPAW()

      else if (dft_program .eq. 'G09') then

        dft_type = 3
        call read_G09()

      else if (dft_program .eq. 'NW') then

        dft_type = 4
        call read_NW()
       
      else if (dft_program .eq. 'VASP') then
        dft_type = 5
        call read_VASP()

      end if

      close(22,status='keep')

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE read_G09()
      use DFT
      use commons, only : GAconv
      use potentials, only: r_ab, namea, nameb, nspec
      implicit none


         read(22,*)nspec
c     Read specifications for DFT runs
         read(22,*)namea
         read(22,*)atn_a
      if (nspec .eq. 2) then
         read(22,*)nameb
         read(22,*)atn_b
      end if

CCCC If needed this can be used for input parameters?
C         read(22,*)sys_DFT(1)   ! 
C         read(22,*)sys_DFT(2)   ! 
C         read(22,*)sys_DFT(3)   ! 
C         read(22,*)sys_DFT(4)   ! 
C         read(22,*)sys_DFT(5)   !
 
         read(22,*)GAconv       ! GAConv
         read(22,*)r_ab         ! Mean bond length

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE read_GPAW()
      USE HECTORSTUFF
      use DFT
      use commons, only : GAconv
      use potentials, only: r_ab, namea, nameb, nspec
      implicit none

         read(22,*)nspec
c     Read specifications for DFT runs
         read(22,*)namea
      if (nspec .eq. 2) then
         read(22,*)nameb
      end if
         read(22,*)scratch      ! Used for the basis set precision
         read(22,*)functional_a ! XC
         read(22,*)functional_b ! Used for the Force convergence criteria on minimisation
         read(22,*)conv_DFT(1)  ! Used for MixerSum (1) beta (2) nmaxold (3) weight
         read(22,*)conv_DFT(2)
         read(22,*)conv_DFT(3)
         read(22,*)sys_DFT(1)   ! Maximum SCF iterations
         read(22,*)sys_DFT(2)   ! Orbital occupation smearing using FermiDirac
         read(22,*)sys_DFT(3)   ! Number of bands occupied
         read(22,*)sys_DFT(4)   ! Eigensolver (See GPAW literature)
         read(22,*)ele_DFT(1)   ! SCF Convergence: Energy
         read(22,*)ele_DFT(2)   ! SCF Convergence: Density
         read(22,*)GAconv       ! GAConv
         read(22,*)r_ab         ! Mean bond length
         read(22,*)openshell    ! Unpaired electrons?
C      if (openshell .eq. 1) then
C         read(22,*)sys_DFT(5)
C      else
C         sys_DFT(5) = '0'
C      end if

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE read_QE()
      use DFT
      use commons, only : GAconv
      use potentials, only: mass_a, mass_b, r_ab, namea, nameb, nspec
      implicit none
      
         read(22,*)nspec
         read(22,*)submit
         read(22,*)namea
         read(22,*)mass_a
      if (nspec .eq. 2) then
         read(22,*)nameb
         read(22,*)mass_b
      end if
         read(22,*)GAconv
         read(22,*)r_ab

      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE read_NW()
      use DFT
      use commons, only : GAconv
      use potentials, only: r_ab, namea, nameb, nspec
      implicit none
 
         read(22,*)nspec
         read(22,*)submit
         read(22,*)namea
      if (nspec .eq. 2) read(22,*) nameb
         read(22,*)GAconv
         read(22,*)r_ab
           
         return
      END  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE read_VASP()
      use DFT
      use commons, only : GAconv
      use potentials, only: mass_a, mass_b, r_ab, namea, nameb, nspec
      implicit none


         read(22,*)nspec
         read(22,*)submit
         read(22,*)namea
         read(22,*)mass_a
      if (nspec .eq. 2) then
         read(22,*)nameb
         read(22,*)mass_b
      end if
         read(22,*)GAconv
         read(22,*)r_ab

      END

