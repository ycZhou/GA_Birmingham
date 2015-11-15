module DFT
! A module to hold all of the DFT optimistion parameters
implicit none
integer :: openshell=0           ! Is a open shell system
integer :: dft_type=0            ! Program to be used for DFT where (1) QE (2) GPAW and (3) G09 (4) NW
integer :: atn_a=0               ! Atomic Number of A
integer :: atn_b=0               ! Atomic Number of B
! This is a local variable to get DFT_program type, then assign integer value
character*(10) dft_program       ! DFT Program to be used.
! All other varialbes are as used by Sven previously
character(len=80) :: functional_a       ! Functional for element a
character(len=80) :: functional_b       ! Functional for element b
!character(len=80) :: namea              ! Element name a
!character(len=80) :: nameb              ! Element name b
character(len=80) :: scratch            ! Specify scratch directory
character(len=80) :: conv_DFT(3)        ! Conv criteria (1) energy (2) force
                                ! (3) optimization step number
character(len=80) :: sys_DFT(9)         ! System variables DFT 
                                 ! (1) E_cut_off (2) occupations
                                 ! (3) smearing (4) degauss
                                 ! (5) number of unpaired electrons
                                 ! (6) total charge (7) vdw switch
                                 ! (8) london_s6 (9) london_rcut
character*(80) :: ele_DFT(3)         ! Electron variables DFT
                                 ! (1) SCF convergence (2) mixing SCF
                                 ! (3) maxstep
!double precision :: r_ab=0d0     ! Radius of bond for initial generation
character(len=80) :: submit      ! Submission line
end module DFT
