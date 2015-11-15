      SUBROUTINE z_heapsort(natoms,cluster,qs)
      implicit none

      integer natoms,ir,l,i,j,k,iarr
      real*8 cluster(3*natoms),qs(natoms),tz,tx,ty,tq

c     Start of Numerical Recipes heap sort

      l=natoms/2+1
      ir=natoms

 10   continue

      if(l.gt.1)then
         l=l-1
         tz=cluster(2*natoms+l)
         ty=cluster(natoms+l)
         tx=cluster(l)
         tq=qs(l)
      else
         tz=cluster(2*natoms+ir)
         ty=cluster(natoms+ir)
         tx=cluster(ir)
         tq=qs(ir)
         cluster(2*natoms+ir)=cluster(2*natoms+1)
         cluster(natoms+ir)=cluster(natoms+1)
         cluster(ir)=cluster(1)
         qs(ir)=qs(1)
         ir=ir-1
         if(ir.eq.1)then
            cluster(2*natoms+1)=tz
            cluster(natoms+1)=ty
            cluster(1)=tx
            qs(1)=tq
            return
         endif
      endif
      i=l
      j=l+l
 20   if(j.le.ir)then
         if(j.lt.ir)then
            if(cluster(2*natoms+j).gt.cluster(2*natoms+j+1))j=j+1
         endif
         if(tz.gt.cluster(2*natoms+j))then
            cluster(2*natoms+i)=cluster(2*natoms+j)
            cluster(natoms+i)=cluster(natoms+j)
            cluster(i)=cluster(j)
            qs(i)=qs(j)
            i=j
            j=j+j	
         else
            j=ir+1
         endif
         go to 20
      endif
      cluster(2*natoms+i)=tz
      cluster(natoms+i)=ty
      cluster(i)=tx
      qs(i)=tq
 
      go to 10

      end
