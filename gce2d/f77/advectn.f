cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advectn (x,x1)

c     ****   compute advection of t, q, etc

      integer nx,nz
      parameter (NX=514,NZ=43)

      real    x(nx,nz),x1(nx,nz)

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg

      real    u1(nx,nz),w1(nx,nz)
      common/bsat1/ u1
      common/badv1/ w1

      real    uu1(nx,nz),ww1(nx,nz)
      common/b2u/ uu1
      common/b2w/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/dumuw/ umd,vmd,wmd

      real    umod(nx,nz),wmod(nx,nz)
      common/bw0/ umod,wmod

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    aconst
      common/bsave/ aconst

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    TQMD(NX,NZ)
c     real    TQMD1(NX,NZ)
      integer i,k
      real    small

      save TQMD
c     save TQMD1

cc    **************************************************************
        do 10 i=1,imax
         dxr(i)=rdx
   10    dxxt(i)=dt*rdx
       do 20 k=2,kles
        dzr(k)=am(k)*rdz/rho(k)
   20   dzzt(k)=dt*am1(k)*rdz
      do 100 k=1,kmax
      do 100 i=1,imax
         umod(i,k)=0.
         wmod(i,k)=0.
        u1(i,k)=.5*(3.*uu1(i,k)-umd(i,k))
        IF (ISMG .NE. 5) THEN
          W1(I,K)=.5*(3.*WW1(I,K)-WMD(I,K))
        ELSE
          W1(I,K)=WW1(I,K)
        ENDIF
  100 continue

          ACONST=0.0
         IF (ISMG .LE. 2) THEN
          CALL TMAXADV (X,SMALL)
           ACONST=ABS(SMALL*1.125)
         ENDIF

      do 150 k=1,kmax
      do 150 i=1,imax
       IF (ISMG .LE. 2) THEN
        X(I,K)=X(I,K)+ACONST
       ENDIF
       x1(i,k)=x(i,k)
  150 continue
      call fadv (x,u1,w1,0)
       call boundy (x,x1)
cc    *********************************************
      call fadvuw (x,x1,u1,w1,umod,wmod)
c      call boundy (umod,wmod)
      call fadv (x,umod,wmod,0)
       call boundy (x,x1)

      IF (ISMG .LE. 2) THEN

       DO 300 K=1,KMAX
       DO 300 I=1,IMAX
        TQMD(I,K)=ACONST
c        TQMD1(I,K)=TQMD(I,K)
c        U1(I,K)=.5*(3.*UU1(I,K)-UMD(I,K))
c        W1(I,K)=.5*(3.*WW1(I,K)-WMD(I,K))
  300 CONTINUE

c      CALL FADV (TQMD,U1,W1,1)

c       do 50 k=2,kles
c        tqmd(1,k)=tqmd(iles,k)
c  50   tqmd(imax,k)=tqmd(2,k)
c      do 60 i=1,imax
c       tqmd(i,1)=tqmd(i,2)
c  60   tqmd(i,kmax)=tqmd(i,kles)

      DO 400 K=1,KMAX
      DO 400 I=1,IMAX
         X(I,K)=X(I,K)-TQMD(I,K)
         X1(I,K)=X1(I,K)-ACONST
  400 CONTINUE
      ENDIF
      return
      end
