      program reopt
      implicit none

      integer mmax,nmax,nspecmax2
      parameter(mmax=20,nmax=200,nspecmax2=9)

      integer natoms,pfunc,idamp,ihard,nspec,nels,m,iprint,i,
     &     nbd(3*nmax),isave(44),iwa(9*nmax),narg,lstr,potlen,olen,
     &     tlen,expo,ifilelen,ilen,dloc,nerr,cutoff

      real*8 off(3*nmax),q(nmax),alpha,r_ab,qa,qb,B_ab,rho_ab,
     &     B_bb,rho_bb,de,re,a2,a3,coeff(11),rcut2,aij(nspecmax2),
     &     vij(nspecmax2),pij(nspecmax2),qij(nspecmax2),r0(nspecmax2),
     &     f,factr,pgtol,g(3*nmax),l(3*nmax),u(3*nmax),dsave(29),
     &     wa(2*mmax*3*nmax+12*nmax+12*mmax*mmax+12*mmax),
     &     no,full,olde,bscale,ram

      character*60 task, csave
      character*40 prog,ifile,ofile,pot,tol
      character*20 inputfile,title
      character*2 el,el1,el2,namea,nameb,elements
      logical      lsave(4),iflag

C ADDED BY ANDREW LOGSDAIL. ANOTHER MESSY FILE
      logical cutoff_gupta
      real*8 cutoff_start(nspecmax2),cutoff_end(nspecmax2)
      real*8 a3a(nspecmax2),a4a(nspecmax2),a5a(nspecmax2)
      real*8 x3x(nspecmax2),x4x(nspecmax2),x5x(nspecmax2)

      external lstr
      
      narg = iargc()

      cutoff_gupta=.false.

c     check number of arguments

      if (narg .lt. 3) then
         call getarg(0,prog)
         ilen = lstr(prog)
         write(*,10)prog(1:ilen)
 10      format('Usage : ',a,' <input file> <potential> <output file>',
     &        ' <tolerance> (optional) cutoff (optional)')
         write(*,15)
 15      format('Default tolerance is 1.d+7',
     &        ', decrease tolerance for more precise optimization')
         stop
      end if
      
      call getarg(1,ifile)
      ilen = lstr(ifile)
 
c     open input file on unit 2

      open(unit=8,status='old',file=ifile(1:ilen),iostat=nerr)

      if (nerr .ne. 0) then
         write(*,'(/,''ERROR, cannot open file : '',a,/)')ifile(1:ilen)
         stop
      end if

      write(*,'(/,''Read input file : '',a,/)')ifile(1:ilen)

c     get potential type

      call getarg(2,pot)
      potlen=lstr(pot)

c     determine potential type

      if (pot(1:potlen) .eq. 'morse') then

         nels = 1
         inputfile = 'morse.in'
         ifilelen = 8
         pfunc = 1
         call read_morse(inputfile,ifilelen,alpha,bscale,el)
         el1 = 'Al'

      else if (pot(1:potlen) .eq. 'mm_C') then
         
         nels = 1
         inputfile = 'mm_C.in'
         ifilelen = 7
         pfunc = 2
         call read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,idamp,
     &        ihard,ram,rcut2,bscale,el,title)
         el1 = 'C'

      else if (pot(1:potlen) .eq. 'mm_Al') then
         
         nels = 1
         inputfile = 'mm_Al.in'
         ifilelen = 8
         pfunc = 2
         call read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,idamp,
     &        ihard,ram,rcut2,bscale,el,title)
         el1 = 'Al'
         
      else if (pot(1:potlen) .eq. 'mm_Yb') then
         
         nels = 1
         inputfile = 'mm_Yb.in'
         ifilelen = 8
         pfunc = 2
         call read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,idamp,
     &        ihard,ram,rcut2,bscale,el,title)
         el1 = 'Yb'

      else if (pot(1:potlen) .eq. 'mm_Pb') then
         
         nels = 1
         inputfile = 'mm_Pb.in'
         ifilelen = 8
         pfunc = 2
         call read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,idamp,
     &        ihard,ram,rcut2,bscale,el,title)
         el1 = 'Pb'

      else if (pot(1:potlen) .eq. 'mm_gen') then
         
         nels = 1
         inputfile = 'mm_gen.in'
         ifilelen = 9
         pfunc = 2
         call read_mmparams(inputfile,ifilelen,de,re,a2,a3,coeff,idamp,
     &        ihard,ram,rcut2,bscale,el,title)

      else if (pot(1:potlen) .eq. 'MgO') then
         
         nels = 2
         inputfile = 'MgO.in'
         ifilelen = 6
         pfunc = 3
         call read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &        B_bb,rho_bb,bscale,namea,nameb)
         el1 = namea
         el2 = nameb

      else if (pot(1:potlen) .eq. 'ZnO') then
         
         nels = 2
         inputfile = 'ZnO.in'
         ifilelen = 6
         pfunc = 3
         call read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &        B_bb,rho_bb,bscale,namea,nameb)
         el1 = namea
         el2 = nameb

      else if (pot(1:potlen) .eq. 'gen_ion') then
         
         nels = 2
         inputfile = 'gen_ion.in'
         ifilelen = 6
         pfunc = 3
         call read_iparams(inputfile,ifilelen,r_ab,qa,qb,B_ab,rho_ab,
     &        B_bb,rho_bb,bscale,namea,nameb)
         el1 = namea
         el2 = nameb

      else if (pot(1:6) .eq. 'gupta_') then
         
         inputfile = pot  
         ifilelen = potlen 
         pfunc = 4
         call read_gupta(inputfile,ifilelen,nspec,aij,vij,pij,qij,r0,
     &        bscale,namea,nameb)
        
         nels = nspec 
         if (nels .eq. 2) then
              el1 = namea
              el2 = nameb
         else
              el1 = namea
         end if
         
      else        
         
         write(*,'("ERROR, Invalid potential type : ",a)')pot(1:potlen)
         stop
         
      end if

c     get output file

      call getarg(3,ofile)
      olen = lstr(ofile)

      open(unit=3,status='unknown',file=ofile(1:olen))

c     Check for Gupta Cutoff

      if (narg .eq. 5) then

         call getarg(5,tol)
         cutoff = lstr(tol)

         if ( tol(1:cutoff) .eq. 'cutoff') then

           call read_cutoff(cutoff_start,cutoff_end,a3a,a4a,a5a,
     &          x3x,x4x,x5x)
         cutoff_gupta = .true.

         end if

      end if

C     Get Tolerance/ Check for Gupta Cutoff

      if (narg .eq. 4) then
         call getarg(4,tol)
         tlen = lstr(tol)

         if ( tol(1:tlen) .eq. 'cutoff') then

           call read_cutoff(cutoff_start,cutoff_end,a3a,a4a,a5a,
     &          x3x,x4x,x5x)
           cutoff_gupta = .true.

         else

C     Cutoff Check done, now just look for tolerance

           dloc = 0

           do i=1,tlen
             if (tol(i:i) .eq. 'd') then
                dloc = i
             end if
           end do

           if (dloc .gt. 0) then

             read(tol(1:dloc-1),'(f4.2)',iostat=nerr)no
             read(tol(dloc+1:tlen),'(i2)',iostat=nerr)expo
 
             if (expo .gt. 0) then
                factr = no*10**(expo-1)
             else
                factr = no/(10**-expo)
             end if

            else

             read(tol(1:tlen),'(f12.3)',iostat=nerr)factr

            end if

           if (nerr .ne. 0) then
              write(*,'("ERROR, invalid tolerance : ",a)')tol(1:tlen)
              stop
           end if

         end if

       else

         factr = 1.d+7

      end if

c     read input file

      if (nels .eq. 1) then

         read(8,'(i3)')natoms
         read(8,*)olde,bscale
          
         do i=1,natoms
            read(8,'(a,1x,f13.8,f13.8,f13.8)')
     &           el,off(i),off(i+natoms),off(i+2*natoms)

        end do
         
      else
         
         if (pfunc .eq. 3) then
         
            read(8,*)natoms
            read(8,*)olde,bscale
            
            do i=1,natoms
               read(8,'(a,1x,f13.8,f13.8,f13.8)') 
     &              el,off(i),off(i+natoms),off(i+2*natoms)
               if (el .eq. el1) then
                  q(i) = qa
               else
                  q(i) = qb
               end if
            end do

         else if (pfunc .eq. 4) then
             
            read(8,*)natoms
            read(8,*)olde,bscale
          
            do i=1,natoms
               read(8,'(a,1x,f13.8,f13.8,f13.8)')
     &              el,off(i),off(i+natoms),off(i+2*natoms)
               if (el .eq. el1) then
                  q(i) = qa
               else
                  q(i) = qb
               end if
            end do

         end if

      endif

      close(8)

c     unscale coordinates
      
      do i=1,natoms
         off(i) = off(i)/bscale
         off(i+natoms) = off(i+natoms)/bscale
         off(i+2*natoms) = off(i+2*natoms)/bscale
      end do 
      
c     perform minimizations
      
c     set up minimisation options
      
      m=5
      pgtol=1.0d-5
      iprint = -1
      
c     no bounds
      
      do i=1,3*natoms
         nbd(i) = 0
      end do
      
      task = 'START'
      
c     start  L-BFGS-B Routine
      
      if (pfunc .eq. 1 ) then
         
c     morse potential
         
         do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)

            if (task(1:2) .eq. 'FG') call morse(natoms,alpha,
     &           off,f,g)
            
         end do
         
      else if ( pfunc .eq. 2 ) then
         
c     murrell-mottram
         
         do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)
            
            if (task(1:2) .eq. 'FG') call murrmott(natoms,idamp,
     &            ihard,de,re,a2,a3,coeff, rcut2,off,f,g)
            
         end do
         
      else if ( pfunc .eq. 3 ) then
         
c     rigid-ion
         
         do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)
            
            if (task(1:2) .eq. 'FG') call ionic(natoms,B_ab,rho_ab,
     &           B_bb,rho_bb,off,f,g,q,qa,qb)
            
         end do
         
      else if ( pfunc .eq. 4 ) then
         
c     gupta potential
         
         do while(task(1:5).eq.'START' .or. task(1:2).eq.'FG'
     &        .or. task(1:5).eq.'NEW_X')

C            write(*,*)task
C            write(*,*)natoms
            
            call setulb(3*natoms,m,off,l,u,nbd,f,g,factr,pgtol,wa,iwa,
     &           task,iprint,csave,lsave,isave,dsave)
            
C            write(*,*)task
C            write(*,*)natoms

            if (task(1:2) .eq. 'FG') call gupta(natoms,nspec,off,q,f,
     &           g,aij,vij,pij,qij,r0,iflag,cutoff_gupta,cutoff_start,
     &           cutoff_end,a3a,a4a,a5a,x3x,x4x,x5x)

         end do
         
      endif
      
c     write new output file - note coordinates are unscaled
      
      write(3,'(i3)')natoms
      write(3,'(f13.6,2x,f4.2)')f,1.00
      
      if (nels .eq. 1) then
      
         do i=1,natoms
            write(3,'(a,1x,3f13.8)')el1,off(i),off(i+natoms),
     &           off(i+2*natoms)
         end do

      else

         if (pfunc .eq. 3) then
         
            do i=1,natoms
               if (q(i) .gt. 0.d0) then
                  write(3,'(a,1x,3f13.8)')el1,off(i),off(i+natoms),
     &           off(i+2*natoms)
               else
                   write(3,'(a,1x,3f13.8)')el2,off(i),off(i+natoms),
     &           off(i+2*natoms)
                end if
             end do

          else

             do i=1,natoms
               if (q(i) .eq. 1.d0) then
                  write(3,'(a,1x,3f13.8)')el1,off(i),off(i+natoms),
     &           off(i+2*natoms)
               else
                   write(3,'(a,1x,3f13.8)')el2,off(i),off(i+natoms),
     &           off(i+2*natoms)
                end if
             end do
             
          end if

       end if

       close(3)

c     write energies to standard output

       write(*,20)olde,f
 20    format('Energy before = ',f13.6,', Energy after = ',f13.6,/)

       stop
       end
                  
      
