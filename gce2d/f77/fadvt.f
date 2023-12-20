


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadvt (x,y,w)

      implicit none
CC    ****   COMPUTE ADVECTION OF TERMINAL VELOCITY

      integer nx,nz
      PARAMETER (NX=514,NZ=43)

      real    x(nx,nz),y(nx,nz),w(nx,nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),
     $   y8(nx),y16(nx),y17(nx),y18(nx),x1(nx),x2(nx),x3(nx),x4(nx),
     $   x5(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,y8,y16,y17,y18,x1,x2,x3,x4,x5

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a0
      save
      DO 100 K=2,KLES
       A0=AM(K)*RDZ*RRHO(K)
      DO 10 I=2,ILES
        X1(I)=-W(I,K)
   10   X2(I)=-W(I,K+1)
      DO 15 I=2,ILES
        Y1(I)=X(I,K+1)
        Y2(I)=X(I,K)
         IF (X2(I) .GE. 0.0) Y1(I)=X(I,K)
         IF (X1(I) .GT. 0.0) Y2(I)=X(I,K-1)
        Y7(I)=-A0*(RHO1(K+1)*X2(I)*Y1(I)-RHO1(K)*X1(I)*Y2(I))
   15  CONTINUE
      DO 100 I=2,ILES
       Y(I,K)=Y(I,K)+Y7(I)
  100 CONTINUE
      RETURN
      END