

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine akcoefRH 
      integer nx,nz
      parameter (NX=514,NZ=43)

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

      real    uuT (nx,nz),wwT(nx,nz),ak1(nx,nz),ak(nx,nz)

      common/b1a/ ak
      common/b2u/ uuT
      common/b2w/ wwT
      common/b2a/ ak1

      real    u1(nx,nz)
      common/badv/ u1

      real    uu1(nx,nz)
      common/bsat/ uu1

      real    ww1(nx,nz)
      common/bsat1/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/dumuw/ umd,vmd,wmd

      real    tbsk(nz),bskt(nz),bsk2t(nz),bsk(nz),bsk4(nz),
     1        bsit(nz),bsi2t(nz),bsi(nz),bsi4(nz)
      common/b4/ tbsk,bskt,bsk2t,bsk,bsk4,bsit,bsi2t,bsi,bsi4

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef22(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef22,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),t(nx),
     $        tp(nx),tm(nx),q(nx),qp(nx),qm(nx),y8(nx),y9(nx),y10(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,t,tp,tm,q,qp,qm,y8,y9,y10

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,km,kp
      real    a111,coef(nz)

      save

c     ********************************
      DO K=1,KMAX
         Y1(K)=DZ/AM(K)
        COEF(K)=.5*CK*CK*Y1(K)*Y1(K)
      enddo
      do 100 k=1,kmax
      do 100 i=1,imax
       IF (IJKADV .EQ. 0) THEN
         UU1(I,K)=UUT(I,K)
         WW1(I,K)=WWT(I,K)
       ELSE
         UU1(I,K)=.5*(3.*UUT(I,K)-UMD(I,K))
         WW1(I,K)=.5*(3.*WWT(I,K)-WMD(I,K))
       ENDIF
  100  ak(i,k)=.5*ak(i,k)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ***   CONVENTIONAL MIXING COEFFICIENT   **************************
c      do 200 k=2,kles
c        kp=k+1
c        km=k-1
c        a111=am(k)*rd2z*rrho(k)*(rho1(kp)-rho1(k))
c        do i=2,iles
c          y1(i)=(ww1(i+1,kp)+ww1(i+1,k)-ww1(i-1,kp)-ww1(i-1,k))*rd4x
c     1         +(uu1(i+1,kp)+uu1(i,kp)-uu1(i+1,km)-uu1(i,km))*rd4z*am(k)
c     2         -a111*(ww1(i,kp)+ww1(i,k))
c          u1(i,k)=coef(k)*y1(i)
c       enddo
c  200  continue
c      do 400 k=2,kles
c       a111=tbsk(k)
c       do i=2,iles
c         ak(i,k)=abs( u1(i,k) )
c          if (ak(i,k) .ge. a111) ak(i,k)=a111
c         ak1(i,k)=ak(i,k)
c       enddo
c  400  continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ***   SCHLESINGER'S MIXING COEFFICIENT   *************************
       A13=1./3.
      do 200 k=2,kles
         kp=k+1
         km=k-1
        a111=a13*am(k)*rd2z*rrho(k)*(rho1(kp)-rho1(k))
        do i=2,iles
           tp(i)=a111*(ww1(i,kp)+ww1(i,k))
          y1(i)=(ww1(i+1,kp)+ww1(i+1,k)-ww1(i-1,kp)-ww1(i-1,k))*rd4x
     1         +(uu1(i+1,kp)+uu1(i,kp)-uu1(i+1,km)-uu1(i,km))*rd4z*am(k)
          y1(i)=y1(i)*y1(i)

          y2(i)=(uu1(i+1,k)-uu1(i,k))*rdx+tp(i)
           y2(i)=y2(i)*y2(i) 
          y3(i)=am(k)*(ww1(i,kp)-ww1(i,k))*rdz+tp(i)
           y3(i)=y3(i)*y3(i)
          u1(i,k)=coef(k)*(2.*(y2(i)+y3(i))+y1(i))
        enddo
  200 CONTINUE

      do 400 k=2,kles
         a111=tbsk(k)
        do i=2,iles
         ak(i,k)=sqrt(u1(i,k))
          if (ak(i,k) .ge. a111) ak(i,k)=a111
         ak1(i,k)=ak(i,k)
        enddo
  400 continue

      call boundy (ak,ak1)

      return
      end
