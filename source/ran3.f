c     Genetic Algorithm Program
c     Develoepd by the Birmingham University Cluster Group, Prof. Roy Johnston
c     Commenting added by Andrew Logsdail, Oct 2010

C     Random number generator
      FUNCTION ran3(idum)
      implicit none
c
c     Numerical Recipes Random Number Generator (ran3)
c

C     Internal integer values
      integer idum
      integer mbig
      integer mseed
      integer mz
      integer i
      integer iff
      integer ii
      integer inext
      integer inextp
      integer k
      integer mj
      integer mk
      integer ma(55)	! Array of size(55)
C     Real values
      real*8 ran3
      real*8 fac
C     Parameters
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
C     Set variables to remain defined after termination of function
      save iff
      save inext
      save inextp
      save ma
C     Define iff value
      data iff /0/

C     Calculate
      if (idum .lt. 0 .or. iff .eq. 0) then
         iff=1
         mj=mseed-iabs(idum)
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 10 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 10      end do
         do 20 k=1,4
            do 30 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 30         end do
 20      end do
         inext=0
         inextp=31
         idum=1
      end if
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac	! Final random output value
     
      return
      end
