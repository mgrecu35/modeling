
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine flxco2(m,n,np,swc,swh,csm,df)

c*****************************************************************

c-----compute the reduction of clear-sky downward solar flux
c     due to co2 absorption.

      implicit none

c-----input parameters

      integer m,n,np
      real    csm(m,n),swc(m,n,np+1),swh(m,n,np+1),cah(22,19)

c-----output (undated) parameter

      real    df(m,n,np+1)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----temporary array

      integer i,j,k,ic,iw,mmm
      real    xx,clog,wlog,dc,dw,x1,x2,y2

c********************************************************************
c-----include co2 look-up table

      include "cah.dat"

      save
 
c********************************************************************
c-----table look-up for the reduction of clear-sky solar
c     radiation due to co2. The fraction 0.0343 is the
c     extraterrestrial solar flux in the co2 bands.
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do k= 2, np+1
       do j= 1, n
        do i= 1, mmm
          xx=1./.3
          clog=log10(swc(i,j,k)*csm(i,j))
          wlog=log10(swh(i,j,k)*csm(i,j))
          ic=int( (clog+3.15)*xx+1.)
          iw=int( (wlog+4.15)*xx+1.)
          if(ic.lt.2)ic=2
          if(iw.lt.2)iw=2
          if(ic.gt.22)ic=22
          if(iw.gt.19)iw=19     
          dc=clog-float(ic-2)*.3+3.
          dw=wlog-float(iw-2)*.3+4.   
          x1=cah(1,iw-1)+(cah(1,iw)-cah(1,iw-1))*xx*dw
          x2=cah(ic-1,iw-1)+(cah(ic-1,iw)-cah(ic-1,iw-1))*xx*dw
          y2=x2+(cah(ic,iw-1)-cah(ic-1,iw-1))*xx*dc
          if (x1.lt.y2) x1=y2
          df(i,j,k)=df(i,j,k)+0.0343*(x1-y2)
        enddo     
       enddo      
      enddo      

      return
      end
