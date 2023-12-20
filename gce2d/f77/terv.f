cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine terv (irsg,rho,fv)
      implicit none
      integer nx,nz,nx10
      parameter (NX=514,NZ=43)
      parameter (nx10=10*nx)
      integer irsg
      real    fv(nz),rho(nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    qrn(nx,nz),qcg(nx,nz),qcs(nx,nz)
      common/b1r/ qrn
      common/b1g/ qcg
      common/b1s/ qcs
      real    ww1(nx,nz)
      COMMON/BTV/ WW1

      real    y1(nx),y2(nx),y3(nx),vr(nx),vs(nx),vg(nx),y4(nx10)
      common/ba/ y1,y2,y3,vr,vs,vg,y4

      real    ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq

      real    zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc
      common/rterv/ zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    vgcr,vscf

      save


c     ******************************************************************
      if(irsg.ne.0) go to 1
      do 10 k=2,kles
      do 10 i=2,iles
       ww1(i,k)=0.
         y1(i)=.5*rho(k)*(qrn(i,k)+qrn(i,k-1))
        if (y1(i) .gt. 1.e-15) then
           vs(i)=sqrt( y1(i) )
           vg(i)=sqrt( vs(i) )
          vr(i)=fv(k)*(vr0+vr1*vg(i)+vr2*vs(i)+vr3*vg(i)*vs(i))
         ww1(i,k)=max(vr(i), 0.E0)
        endif
   10 continue
      return
    1 if(irsg.ne.1) go to 2
      do 20 k=2,kles
       vscf=vsc*fv(k)
      do 20 i=2,iles
       ww1(i,k)=0.
         y1(i)=.5*rho(k)*(qcs(i,k)+qcs(i,k-1))
        if (y1(i) .gt. 1.e-15) then
         vs(i)=vscf*y1(i)**bsq
         ww1(i,k)=max(vs(i), 0.E0)
        endif
   20 continue
      return
    2 do 30 k=2,kles
       vgcr=vgc*fv(k)
      do 30 i=2,iles
       ww1(i,k)=0.
         y1(i)=.5*rho(k)*(qcg(i,k)+qcg(i,k-1))
        if (y1(i) .gt. 1.e-15) then
          vg(i)=vgcr*y1(i)**bgq
         ww1(i,k)=max(vg(i), 0.E0)
        endif
   30 continue
      return
      end
