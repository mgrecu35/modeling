cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advectak (x,x1)

      implicit none
CC    ****   COMPUTE ADVECTION OF K COEFFICIENT
      integer nx,nz
      PARAMETER (NX=514,NZ=43)

      real    x(nx,nz),x1(nx,nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    u1(nx,nz),w1(nx,nz)
      common/bsat/ u1
      common/bsat1/ w1

      real    uu1(nx,nz),ww1(nx,nz)
      common/b2u/ uu1
      common/b2w/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz),umod(nx,nz),wmod(nx,nz)
      common/dumuw/ umd,vmd,wmd
      common/bw0/ umod,wmod

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k

      save


CC    **************************************************************
        DO 10 I=1,IMAX
         DXR(I)=RDX
   10    DXXT(I)=DT*RDX
       DO 20 K=2,KLES
        DZR(K)=AM(K)*RDZ/RHO(K)
   20   DZZT(K)=DT*AM1(K)*RDZ
      DO 100 K=1,KMAX
      DO 100 I=1,IMAX
        UMOD(I,K)=0.
        WMOD(I,K)=0.
       U1(I,K)=.5*(3.*UU1(I,K)-UMD(I,K))
       W1(I,K)=.5*(3.*WW1(I,K)-WMD(I,K))
  100  X1(I,K)=X(I,K)
      DO 150 K=1,KMAX
      DO 150 I=1,IMAX
        IF (X(I,K) .LT. 1.E-5) X(I,K)=0.0
       X1(I,K)=X(I,K)
  150 CONTINUE
      CALL FADVAK (X,U1,W1)
       CALL BOUNDY (X,X1)
       CALL FADVUW (X,X1,U1,W1,UMOD,WMOD)
       CALL FADVAK (X,UMOD,WMOD)
       CALL BOUNDY (X,X1)
      DO 200 K=1,KMAX
      DO 200 I=1,IMAX
        IF (X(I,K) .LT. 1.E-5) X(I,K)=0.
  200 CONTINUE
      RETURN
      END
