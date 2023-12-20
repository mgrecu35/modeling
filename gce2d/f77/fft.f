cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fft (iflg,l)

      implicit none

      integer nx,nz,nxi,nxm,nzm,nxq,nb,iflg,l
      parameter (NX=514,NZ=43)
      parameter (nxi=nx-2,nxm=nx-1,nzm=nz-1,nxq=nxi/2+1,nb=10*nx)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nb)
      common/ba/ y1,y2,y3,y4,y5,y6,y7

      real    ar(nxm,nzm),ai(nxm,nzm)
      common/rfft/ ar,ai

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer iv(nxi)
      real    s(nxq)
      integer i,i1,i2,id1,ie,ivi,ix,j,jj,jmm,jmmn,k,ks
      integer limit,lo2,lo4,mm,n2
      real    cosv,pint,sinv

      save

c
      if(iflg.eq.1) go to 100
      ie=0
      i=1
    6 ie=ie+1
      i=i+i
      if(i.lt.l) go to 6
      pint=8.*atan(1.)/float(l)
      limit=l/4+1
      do 7 i=1,limit
    7  s(i)=sin((i-1)*pint)
      iv(1)=0
      do 5 i=2,l
       ivi=0
       i1=i-1
       i2=i1/2
       do 8 j=1,ie
        ivi=ivi+ivi+i1-i2-i2
        i1=i2
        i2=i2/2
    8  continue
       iv(i)=ivi
    5 continue
      return
  100  lo2=l/2
       lo4=l/4
       n2=lo2
      do 1500 i=1,n2
       do 150 k=1,kl2
        y1(k)=ar(i,k)+ar(i+n2,k)
        y2(k)=ai(i,k)+ai(i+n2,k)
        y3(k)=ai(i,k)-ai(i+n2,k)
        y4(k)=ar(i,k)-ar(i+n2,k)
        ar(i,k)=y1(k)
        ai(i,k)=y2(k)
        ai(i+n2,k)=y3(k)
        ar(i+n2,k)=y4(k)
  150 continue
 1500 continue
       mm=l
       jj=1
  500  mm=mm/2
       n2=n2/2
       jj=jj+jj
       ks=-1
      do 3000 j=1,jj
       ks=ks+1
       ix=iv(ks+ks+1)
       jmm=mm*(j-1)
       jmmn=jmm+n2
       if(lo4-ix) 20,15,15
   15   sinv=s(ix+1)
        cosv=s(lo4-ix+1)
       go to 25
   20   sinv=s(lo2-ix+1)
        cosv=-s(ix-lo4+1)
   25  continue
       do 300 i=1,n2
        do 30 k=1,kl2
         y1(k)=ar(i+jmmn,k)*cosv-ai(i+jmmn,k)*sinv
         y2(k)=ar(i+jmmn,k)*sinv+ai(i+jmmn,k)*cosv
         y3(k)=ar(i+jmm,k)+y1(k)
         y4(k)=ai(i+jmm,k)+y2(k)
         y5(k)=ar(i+jmm,k)
         y6(k)=ai(i+jmm,k)
         ar(i+jmmn,k)=y5(k)-y1(k)
         ai(i+jmmn,k)=y6(k)-y2(k)
         ar(i+jmm,k)=y3(k)
   30    ai(i+jmm,k)=y4(k)
  300  continue
 3000 continue
      if(n2.ge.2) go to 500
      do 600 i=1,l
        id1=iv(i)+1
        if(i.ge.id1) go to 600
       do 60 k=1,kl2
         y1(k)=ar(id1,k)
         y2(k)=ai(id1,k)
         ar(id1,k)=ar(i,k)
         ai(id1,k)=ai(i,k)
         ar(i,k)=y1(k)
         ai(i,k)=y2(k)
   60  continue
  600 continue
      return
      end
