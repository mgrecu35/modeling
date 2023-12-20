Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmpf
      implicit none
c     ******   compute force terms   ***********************************
c     ******   update w, u and v without pressure force

      integer nx,nz,nx10
      parameter (NX=514,NZ=43)
      parameter (nx10=10*nx)

      integer LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    uu1(nx,nz),vv1(nx,nz),ww1(nx,nz),
     $        u  (nx,nz),v  (nx,nz),w  (nx,nz)
      common/b1u/ u
      common/b1v/ v
      common/b1w/ w
      common/b2u/ uu1
      common/b2v/ vv1
      common/b2w/ ww1

      real    f(nx,nz),w1(nx,nz)
      common/bsat/ f
      common/bsat1/ w1

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx10)
      common/ba/ y1,y2,y3,y4,y5,y6,y7

      real UTEM_DT(NX,NZ),VTEM_DT(NX,NZ)
      COMMON/DEBUG_UUVV/ UTEM_DT,VTEM_DT

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a1,a2,rd2t

      save

c     ********************************
      rd2t=1./d2t
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO 10 K=2,KLES
        Y1(K)=AM(K)*RDZ*RRHO(K)
        IF (LIPPS .eq. 0) THEN
          Y2(K)=DZ2/(CP*TB(K))
        else
          y2(k)=dz2/cp
        endif
   10   F(IMAX,K)=F(2,K)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 1000 k=2,kles
        do 20 i=2,imax
   20    y3(i)=f(i,k)
       do 100 i=2,iles
       f(i,k)=((y3(i+1)-y3(i))*rdx
     1        +(w1(i,k+1)*rho1(k+1)-w1(i,k)*rho1(k))*y1(k)
     2       +((uu1(i+1,k)-uu1(i,k))*rdx
     3       +(ww1(i,k+1)*rho1(k+1)-ww1(i,k)*rho1(k))*y1(k))*rd2t)*y2(k)

          UTEM_DT(I,K)=UU1(I,K)

         y4(i)=u(i,k)+eps*(-2.*u(i,k)+uu1(i,k))
        u(i,k)=uu1(i,k)+y3(i)*d2t
        uu1(i,k)=y4(i)
c
         y5(i)=w(i,k)+eps*(-2.*w(i,k)+ww1(i,k))
        w(i,k)=ww1(i,k)+w1(i,k)*d2t
  100   ww1(i,k)=y5(i)
 1000 continue
c   ****   set boundary conditions  ************************************
      IF (LIPPS .eq. 0) THEN
        a1=dz*ba(2)/(cp*tb(2)*am1(2))
        a2=dz*bb(kles)/(cp*tb(kles)*am1(kles))
      ELSE
        A1=DZ*BA(2)/(CP*AM1(2))
        A2=DZ*BB(KLES)/(CP*AM1(KLES))
      ENDIF
       do 60 i=2,iles
        f(i,2)=f(i,2)+w1(i,2)*a1
   60   f(i,kles)=f(i,kles)-w1(i,kmax)*a2
C
      return
      end
