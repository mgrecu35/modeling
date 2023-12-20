cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ttop (twc,ktop,kbase)

      implicit none
      integer nz,nzp,nzm
      parameter (NZ=43,nzp=nz+1,nzm=nz-1)

      integer kbase,ktop
      real    twc(nz)

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer k,k1,ktemp

      save


c
        ktop=2
        do 10 k=2,nzm
         if (twc(k).lt.1.e-5) goto 10
          ktemp=k
         if(ktop.lt.ktemp) ktop=ktemp
   10   continue
         kbase=nzm
        do 20 k1=2,nzm
         k=nzp-k1
         if (twc(k).lt.1.e-5) goto 20
          ktemp=k
          if (kbase.gt.ktemp) kbase=ktemp
   20    continue
      if(ktop.lt.kbase)then
        ktop=1
        kbase=2
      endif
      return
      end   
