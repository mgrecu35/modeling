
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine domain (x,i1,i2,k1,k2)

      implicit none

      integer nx,nz
      parameter (NX=514,NZ=43)

      integer i1,i2,k1,k2
      real    x(nx,nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    y

      save

      i1=0
      y=0.
      do 100 i=2,iles
       do 10 k=2,kles
   10   y=max(y,x(i,k))
       if(y .gt. 0.) then
        i1=i
        go to 15
       endif
  100 continue
   15 if(i1 .eq. 0) then
       i1=2
       i2=2
       k1=2
       k2=2
       return
      endif
      y=0.
      do 200 i=iles,2,-1
       do 20 k=2,kles
   20  y=max(y,x(i,k))
       if(y .gt. 0.) then
        i2=i
        go to 25
       endif
  200 continue
   25 y=0.
      do 300 k=2,kles
       do 30 i=2,iles
   30  y=max(y,x(i,k))
       if(y .gt. 0.) then
        k1=k
        go to 35
       endif
  300 continue
   35 y=0.
      do 400 k=kles,2,-1
       do 40 i=2,iles
   40  y=max(y,x(i,k))
       if(y .gt. 0.) then
        k2=k
        go to 45
       endif
  400 continue
   45 i1=max(2,i1-3)
      i2=min(iles,i2+3)
      k1=max(2,k1-3)
      k2=min(kles,k2+3)
      return
      end
