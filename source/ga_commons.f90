module commons
implicit none
! Parameters from ga.param
! Maximum dimensions of various arrays.
! Probably better to use dynamic allocation.
integer :: nmax=320
integer :: ncmax=100
integer :: maxgen=500
integer :: fnamelen=20
integer :: nspecmax2=9
integer :: maxtoklen=20
integer :: maxtokens=40
integer :: tokbufflen=800!maxtoklen*maxtokens

!     integers

integer natoms
integer :: na = 0
integer :: nb = 0
integer nclust
integer ngen
integer noff
integer pfunc
integer idum
integer flen
integer mutant(100) !ncmax
integer nstride
integer meswaps		! Used for mutate_exchange, defines number of swaps to be performed
integer processors        ! Currently will be read, but unused in the ga due to OMP problems
integer tcount ! Terminate if no improvement for tcount generations
integer tsize  ! Tournament Size
integer ptype     ! Potential Type
integer mtype     ! Mating Type
integer stype     ! Selection Type
integer ftype     ! Fitness Type
integer muttype   ! Mutation Type
integer nels      ! Number of Elements
integer mscheme   ! Mutation Scheme

!     reals

real*8 pi
real*8 range
real*8 mrate
real*8 merate		! Frequency with which to perform mutate exchange on top of other mutations
real*8 gmin
real*8 epred! Energy threshold for predator
double precision :: GAconv=0d0   ! GA convergence criteria

!     character

character*(20) fname !fnamelen

!     logical
 
logical tflag
logical highemut ! Use High Energy Mutants
logical remove ! Remove similar clusters	
logical write_clust       ! Save clusters
logical write_energy      ! Save energies
logical write_stats       ! Save statistics
logical pred ! High energy predator
logical cutoff_gupta      ! Use gupta potential cutoff
logical restart	! Using restart

!     array declarations
! Allocate these dynamically one day.
real*8 pop(96100)!((3*nmax+1)*ncmax)
real*8 energy(100)!(ncmax)
real*8 fit(100)!(ncmax)
real*8 q(32000)!(nmax*ncmax)
real*8 qp(640)!(2*nmax)
real*8 pair(1922)!(6*nmax+2)
real*8 off(961)!(3*nmax+1)
real*8 qoff(320)!(nmax)

end module commons

module potentials
implicit none
!     potential parameters

!     common

real*8 bscale
character*20 inputfile
character*2 el
!     Morse

real*8 alpha
      
!     Murrell-Mottram

integer idamp
integer ihard
real*8 de
real*8 re
real*8 a2
real*8 a3
real*8 coeff(11)
real*8 ram
real*8 rcut2
character*6 title

!     Rigid Ion

double precision :: r_ab=0d0     ! Radius of bond for initial generation
real*8 qa
real*8 qb
real*8 B_ab
real*8 rho_ab
real*8 B_bb
real*8 rho_bb
character(len=2) :: namea              ! Element name a
character(len=2) :: nameb              ! Element name b

!     Gupta 

integer nspec
!MTO
!Originally these were all icompiled length nspecmax2.
!For now, I'll make them all fixed-length.
real*8 aij(9)
real*8 vij(9)
real*8 pij(9)
real*8 qij(9)
real*8 r0(9)

!     Silica
real*8 ru_sisi
real*8 ru_sio
real*8 ru_oo
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
double precision :: mass_a=0d0   ! Mass of element a
double precision :: mass_b=0d0   ! Mass of element b
end module potentials
