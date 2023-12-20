cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advwuv (nudge)
      integer nx,nz,nx10,itt
      parameter (NX=514,NZ=43,ITT=244,ITT3=ITT*3)
      parameter (nx10=10*nx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      COMMON/SFLUXS/ SUW,SVW,SWT,SWQ
      integer itoga,ISFC,ICE,ICE2
      common/itoga/ itoga,ISFC,ICE,ICE2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/cpsbm/ icps_bm,iexplicit

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      integer irf,iadvh,irfg,idq,ismg
      common/bcor/ irf,iadvh,irfg,idq,ismg

      real    fcor
      common/bcor1/ fcor

      common/add_f/ i_four

      integer n,isec,nran
      common/bstart/ n,isec,nran

      real    dpt(nx,nz),dqv(nx,nz),qcl(nx,nz),qrn(nx,nz),
     $        qci(nx,nz),qcs(nx,nz),qcg(nx,nz),
     $        u  (nx,nz),v  (nx,nz),w  (nx,nz),ak (nx,nz),
     $        uu1(nx,nz),vv1(nx,nz),ww1(nx,nz)
      common/b1t/ dpt
      common/b1q/ dqv
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1i/ qci
      common/b1s/ qcs
      common/b1g/ qcg
      common/b1u/ u
      common/b1v/ v
      common/b1w/ w
      common/b1a/ ak
      common/b2u/ uu1
      common/b2v/ vv1
      common/b2w/ ww1

      real    u1(nx,nz)
      common/bsat/ u1

      real    w1(nx,nz)
      common/bsat1/ w1

      real    v1(nx,nz)
      common/badv/ v1

      real    a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
      common/o4/ a4k,a4h,d58x,d16x,d48x,d516x,d32x,d96x

      real    tbsk(nz),bskt(nz),bskt2(nz),bsk(nz),bsk4(nz),
     1        bsit(nz),bsit2(nz),bsi(nz),bsi4(nz)
      common/b4/ tbsk,bskt,bskt2,bsk,bsk4,bsit,bsit2,bsi,bsi4

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/bb6/ tls1,tls2,qls1,qls2,tls3,tls4,qls3,qls4,sft,sfq,wbt,
     $            wb_6h,ub_6h,ubt,q1_6h,q1t,q2_6h,q2t,vb_6h,vbt
      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx10)
      common/ba/ y1,y2,y3,y4,y5,y6,y7
      real    y11(nx),y22(nx),y33(nx),y44(nx),y55(nx),y66(nx)
      real    y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/bls/y0,ths,qs,ts,pss
      real    rfa(nz),rfa1(nz),cntc(nz),cgwd
      common/damp/ rfa,rfa1,cntc,cgwd

      common/bls3/ factor_nuding
      common/bls4/ ubi(nz),vbi(nz),ub_2h(nz,itt3),vb_2h(nz,itt3)
C
CC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFF_VX(NX,NZ) :  V-TURBULENT IN X-COMPNT                        C
C     DIFF_VZ(NX,NZ) :  V-TURBULENT IN Z-COMPNT                        C
C     GRID_VX(NX,NZ) :  V-GRID SCALE TRANSPT IN X-COMPNT               C
C     GRID_VZ(NX,NZ) :  V-GRID SCALE TRANSPT IN Z-COMPNT               C
C     DIFF_NV(NX,NZ) :  V-NUMERICAL FILTER                             C
C     V_LARGE(NX,NZ) :  LARGE-SCALE FORCING IN V WIND                  C
C     PRE_V(NX,NZ)   :  HORIZONTAL PREES FORCING                       C
C     DT_VWIND(NX,NZ):  LOCAL TIME CHANGE TERM                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFF_UX(NX,NZ) :  U-TURBULENT IN X-COMPNT                        C
C     DIFF_UZ(NX,NZ) :  U-TURBULENT IN Z-COMPNT                        C
C     GRID_UX(NX,NZ) :  U-GRID SCALE TRANSPT IN X-COMPNT               C
C     GRID_UZ(NX,NZ) :  U-GRID SCALE TRANSPT IN Z-COMPNT               C
C     DIFF_NU(NX,NZ) :  U-NUMERICAL FILTER                             C
C     U_LARGE(NX,NZ) :  LARGE-SCALE FORCING IN U WIND                  C
C     PRE_U(NX,NZ)   :  HORIZONTAL PREES FORCING                       C
C     DT_UWIND(NX,NZ):  LOCAL TIME CHANGE TERM                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real            DIFF_VX(NX,NZ),DIFF_VZ(NX,NZ),GRID_VX(NX,NZ),
     1                GRID_VZ(NX,NZ),DIFF_NV(NX,NZ),V_LARGE(NX,NZ),
     2                PRE_V(NX,NZ),DT_VWIND(NX,NZ)
      COMMON/DEBUG_V/ DIFF_VX,DIFF_VZ,GRID_VX,
     1                GRID_VZ,DIFF_NV,V_LARGE,
     2                PRE_V,DT_VWIND

      real            DIFF_UX(NX,NZ),DIFF_UZ(NX,NZ),GRID_UX(NX,NZ),
     1                GRID_UZ(NX,NZ),DIFF_NU(NX,NZ),U_LARGE(NX,NZ),
     2                PRE_U(NX,NZ),DT_UWIND(NX,NZ)
      COMMON/DEBUG_U/ DIFF_UX,DIFF_UZ,GRID_UX,
     1                GRID_UZ,DIFF_NU,U_LARGE,
     2                PRE_U,DT_UWIND

      real UTEM_DT(NX,NZ),VTEM_DT(NX,NZ)
      COMMON/DEBUG_UUVV/ UTEM_DT,VTEM_DT

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,km,kmm,kp,im
      real    a0000,c33k,rdxz4,ub2,BARFCT,rd2t,ub12,vb12

      save

c     ******   compute w1   ********************************************
      A0000=1.
       IF (IMLIFTING .EQ. 0) A0000=0.
       a_irf=0.
         if (irf .eq. 1) a_irf=1.
c
      do 10 i=2,iles
       y1(i)=0.
       y3(i)=0.
   10 continue
      do 1000 k=2,kles
        kp=k+1
c        a11k=AM(K)*RDZ2
c        a22k=am1(k)*rrho1(k)
       do 100 i=2,iles
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDANCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        y2(i)=rho(k)*(w(i,kp)+w(i,k))**2*rd4z
        y4(i)=rho(k)*(-am(k)*
     1        (2.*ak(i,k)+bsk(k))
c    1        (2.*ak(i,k)+.5*(bsk(k)+bsk(kp)))
     2        *(ww1(i,kp)-ww1(i,k))*rdz2
     3        +c3(k)*ak(i,k)**2*rdz)
       w1(i,k)=am1(k)*(y1(i)-y2(i)+y3(i)-y4(i))*rrho1(k)
        y1(i)=y2(i)
        y3(i)=y4(i)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  100  continue
 1000 continue
c      do 115 i=2,iles
c        w1(i,kmax)=am1(kmax)*(y2(i)+y4(i))*rrho1(kmax)
c  115 continue

cccccc  tao (11-13-97)
        if (ijkadv .eq. 0 .and. I_FOUR .EQ. 1) then
cc        IF (I_FOUR .EQ. 1) THEN
c
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c
c          do i=2,iles
c             do k=3,kl2
c               y1(k)=(ww1(i,k+2)-3.*(ww1(i,k+1)-ww1(i,k))-ww1(i,k-1))
c             enddo
c             do k=4,kl2
c               w1(i,k)=w1(i,k)+r2dt2*(y1(k-1)-y1(k))
c             enddo
c             w1(i,3)=w1(i,3)+r2dt8*(ww1(i,4)-2.*ww1(i,3))
c             w1(i,kles)=w1(i,kles)+r2dt8*(-2.*ww1(i,kles)+ww1(i,kl2))
c          enddo
c
c          IF (I_FOUR .EQ. 1) THEN
c
               r2dt2=0.005/dt
c
          do i=2,iles
             do k=2,kles
               y1(k)=ww1(i,k+1)-2.*ww1(i,k)+ww1(i,k-1)
             enddo
             y1(2)=y1(3)
             y1(kmax)=y1(kles)
             do k=3,kles
               w1(i,k)=w1(i,k)-r2dt2*(y1(k+1)-2.*y1(k)+y1(k-1))
             enddo
          enddo

c
c
c        ELSE
c
c          y1(2)=0.
c          y1(kmax)=0.
c               r2dt8=0.008/d2t
c          do i=2,iles
c
c            do k=3,kles
c               y1(k)=ww1(i,k)-ww1(i,k-1)
c            enddo
c            do k=3,kles
c               w1(i,k)=w1(i,k)-r2dt8*(y1(k)-y1(k+1))
c            enddo

c          enddo
C
C

cc        ENDIF
        endif
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       IF (IMLIFTING .EQ. 1) THEN
        do 30 k=2,kmax
   30    y1(k)=0.
        do 40 i=2,iles
        do 40 k=2,kmax
   40    y1(k)=y1(k)+w1(i,k)
        do 45 k=2,kmax
   45     y1(k)=y1(k)*ril2
        do 50 k=2,kmax
        do 50 i=2,iles
   50    w1(i,k)=w1(i,k)-y1(k)
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        rdxz4=rd4z*rdx
      do 1500 k=3,kles
        kp=k+1
        km=k-1
        kmm=k-2
c        a11k=AM1(K)*RD2Z
c        a22k=am1(k)*rdxz4
       if (iadvh .eq. 2) then
       do 120 i=2,imax
         y6(i)=ak(i,k)+ak(i,km)+ak(i-1,k)+ak(i-1,km)
        y1(i)=(u(i,k)+u(i,km))*(w(i,k)+w(i-1,k))*rd4x
c    1        -(y6(i)+bsi4(k))*(ww1(i,k)-ww1(i-1,k))*r4dx2
     1        -y6(i)*(ww1(i,k)-ww1(i-1,k))*r4dx2
     2        -am1(k)*y6(i)*(uu1(i,k)-uu1(i,km))*rd4z*rdx
  120  continue
       do 130 i=2,iles
  130   w1(i,k)=w1(i,k)+y1(i)-y1(i+1)
c       y1(1)=(ww1(3,k)-3.*(ww1(2,k)-ww1(1,k))-ww1(il2,k))
c       do 125 i=2,il2
c  125   y1(i)=(ww1(i+2,k)-3.*(ww1(i+1,k)-ww1(i,k))-ww1(i-1,k))
c       do 135 i=2,il2
c  135   w1(i,k)=w1(i,k)+.5*(a4k(k)+a4k(km))*(y1(i-1)-y1(i))
c       w1(iles,k)=w1(iles,k)+.5*(a4k(k)+a4k(km))*(y1(il2)-y1(1))
       y1(1)=a4k(k)*(ww1(3,k)-3.*(ww1(2,k)-ww1(1,k))-ww1(il2,k))
       do 125 i=2,il2
  125   y1(i)=a4k(k)*(ww1(i+2,k)-3.*(ww1(i+1,k)-ww1(i,k))-ww1(i-1,k))
       do 135 i=2,il2
  135   w1(i,k)=w1(i,k)+y1(i-1)-y1(i)
       w1(iles,k)=w1(iles,k)+y1(il2)-y1(1)
       else
c     ***   4-th order horizontal advection
       y3(1)=-d96x*(u(1,kmm)+u(1,kp))*(w(il2,k)+w(1,k))
       do 140 i=2,imax
         y6(i)=ak(i,k)+ak(i,km)+ak(i-1,k)+ak(i-1,km)
        y1(i)=-y6(i)*((ww1(i,k)-ww1(i-1,k))*r4dx2
     1        +am1(k)*(uu1(i,k)-uu1(i,km))*rd4z*rdx)
        y3(i)=-d96x*(u(i,kmm)+u(i,kp))*(w(i-1,k)+w(i,k))
  140  continue
       do 150 i=2,iles
  150   w1(i,k)=w1(i,k)+y1(i)-y1(i+1)
       y1(1)=(u(2,km)+u(2,k))*(d516x*(w(1,k)+w(2,k))
     1                         -d32x*(w(il2,k)+w(3,k)))
     2       +a4k(k)*
c    2       +.5*(a4k(k)+a4k(km))*
     3               (ww1(3,k)-3.*(ww1(2,k)-ww1(1,k))-ww1(il2,k))
       do 160 i=2,il2
        y1(i)=(u(i+1,km)+u(i+1,k))*(d516x*(w(i,k)+w(i+1,k))
     1                              -d32x*(w(i-1,k)+w(i+2,k)))
     2        +a4k(k)*
c    2        +.5*(a4k(k)+a4k(km))*
     3             (ww1(i+2,k)-3.*(ww1(i+1,k)-ww1(i,k))-ww1(i-1,k))
  160  continue
       do 170 i=2,il2
  170   w1(i,k)=w1(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)
       w1(iles,k)=w1(iles,k)+y1(il2)-y1(1)+y3(il2)-y3(3)
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO 175 I=2,ILES
  175   W1(I,K)=W1(I,K)
c    1         -a000*a11k*Y0(I)*(WB(KP)*W(I,KP)-WB(KM)*W(I,KM))
     1         -a0000*am1(k)*y0(i)*(wb(kp)*w(i,kp)-wb(km)*w(i,km))*rd2z
 1500 continue
      do i=2,iles
         w1(i,2)=0.
         w1(i,kmax)=0.
      enddo
c     ******   compute buoyancy
      do 180 i=2,iles
        y1(i)=490.*(dpt(i,1)/tb(1)+.61*dqv(i,1)-qcl(i,1)
     1        -qrn(i,1)-qci(i,1)-qcs(i,1)-qcg(i,1)+a0000*sq(1))
  180  continue
      do 1800 k=2,kmax
       do 190 i=2,iles
         y2(i)=490.*(dpt(i,k)/tb(k)+.61*dqv(i,k)-qcl(i,k)
     1         -qrn(i,k)-qci(i,k)-qcs(i,k)-qcg(i,k)+a0000*sq(k))
        w1(i,k)=w1(i,k)+y1(i)+y2(i)
        w1(i,k)=w1(i,k)-a_irf*rfa1(k)*ww1(i,k)
  190    y1(i)=y2(i)
 1800 continue
c     ******   compute u1   *******************************************
      do i=2,iles
         y1(i)=0.
         y2(i)=0.
         y3(i)=0.
         y4(i)=0.
         y5(i)=0.
         y6(i)=0.
         y7(i)=0.
         y11(i)=0.
         y22(i)=0.
         y33(i)=0.
         y44(i)=0.
         y55(i)=0.
         y66(i)=0.
      enddo
C
C-----put the surface momentum fluxes in here
C
      IF (ITOGA .EQ. 1 .AND. ISFC .EQ. 1) THEN
        BARFCT=.5*RHO1(2)*100.*100.*RDZ  !convert from mks to cgs
        DO I=2,ILES
            im=i-1
            if (i.eq.2) im=iles
          y3(i)=(suw(i)+suw(im))*barfct
          y33(i)=(svw(i)+svw(im))*barfct
        ENDDO
      ENDIF
       rdxz4=rd4x*rdz
      do 2000 k=3,kles
       km=k-1
c        a11_4z2=am1(k)*r4dz2
c        a11_rhoz=rho1(k)*rd4z
c        a11z=am(km)*rrho(km)
c        BASE_DIFF=rho1(k)*am1(k)*bsk4(k)*r4dz2
c
        ub12=ub(k)+ub(km)
        vb12=vb(k)+vb(km)
c
       do 200 i=2,iles
         y7(i)=ak(i,k)+ak(i,km)+ak(i-1,k)+ak(i-1,km)
c
c         y2(i)=rho1(k)*((w(i,k)+w(i-1,k))*(u(i,k)+u(i,km))*rd4z
c     1        -am1(k)*r4dz2*(y7(i)*(uu1(i,k)-uu1(i,km))
c     2        +bsk4(k)*(uu1(i,k)-uu1(i,km)-ub1(k)+ub1(km)))
c     3        -y7(i)*(ww1(i,k)-ww1(i-1,k))*rdxz4)
c
c         y22(i)=rho1(k)*((w(i,k)+w(i-1,k))*(v(i,k)+v(i,km))*rd4z
c     1        -am1(k)*r4dz2*(y7(i)*(vv1(i,k)-vv1(i,km))
c     2        +bsk4(k)*(vv1(i,k)-vv1(i,km)-vb1(k)+vb1(km))))
c

         y2(i)=rho1(k)*(w(i,k)+w(i-1,k))*(u(i,k)+u(i,km))*rd4z
         Y6(I)=RHO1(K)*(W(I,K)+W(I-1,K))*(U(I,K)+U(I,KM)-UB12)*RD4Z

c         y4(i)=-rho1(k)*(am1(k)*r4dz2*y7(i)*(uu1(i,k)-uu1(i,km))
c     1                               +y7(i)*(ww1(i,k)-ww1(i-1,k))*rdxz4)
c         y6(i)=-rho1(k)*am1(k)*r4dz2*bsk4(k)*(uu1(i,k)-uu1(i,km)
c     1                                       -ub1(k)+ub1(km))

         y4(i)=rho1(k)*(-am1(k)*r4dz2*(y7(i)*(uu1(i,k)-uu1(i,km))
     1                  +bsk4(k)*(uu1(i,k)-uu1(i,km)-ub1(k)+ub1(km)))
     2                  -y7(i)*(ww1(i,k)-ww1(i-1,k))*rdxz4)

c

         y22(i)=rho1(k)*(w(i,k)+w(i-1,k))*(v(i,k)+v(i,km))*rd4z
         Y66(I)=RHO1(K)*(W(I,K)+W(I-1,K))*(V(I,K)+V(I,KM)-VB12)*RD4Z

c         y44(i)=-rho1(k)*am1(k)*r4dz2*y7(i)*(vv1(i,k)-vv1(i,km))
c         y66(i)=-rho1(k)*am1(k)*r4dz2*bsk4(k)*(vv1(i,k)-vv1(i,km)
c     1                                  -vb1(k)+vb1(km))

         y44(i)=rho1(k)*(-am1(k)*r4dz2*(y7(i)*(vv1(i,k)-vv1(i,km))
     1                   +bsk4(k)*(vv1(i,k)-vv1(i,km)-vb1(k)+vb1(km))))

c
        u1(i,km)=am(km)*rrho(km)*(y1(i)-y2(i)+y3(i)-y4(i))
        v1(i,km)=am(km)*rrho(km)*(y11(i)-y22(i)+y33(i)-y44(i))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         GRID_UZ(I,KM)=GRID_UZ(I,KM)+am(km)*rrho(km)*(Y1(I)-Y2(I))
c         DIFF_NU(I,KM)=DIFF_NU(I,KM)+am(km)*rrho(km)*(Y3(I)-Y4(I))
c         DIFF_UZ(I,KM)=DIFF_UZ(I,KM)+am(km)*rrho(km)*(Y5(I)-Y6(I))
C
c         GRID_VZ(I,KM)=GRID_VZ(I,KM)+am(km)*rrho(km)*(Y11(I)-Y22(I))
c         DIFF_NV(I,KM)=DIFF_NV(I,KM)+am(km)*rrho(km)*(Y33(I)-Y44(I))
c         DIFF_VZ(I,KM)=DIFF_VZ(I,KM)+am(km)*rrho(km)*(Y55(I)-Y66(I))
C
         GRID_UZ(I,KM)=am(km)*rrho(km)*(Y1(I)-Y2(I))
         DIFF_NU(I,KM)=am(km)*rrho(km)*(Y3(I)-Y4(I))
         DIFF_UZ(I,KM)=am(km)*rrho(km)*(Y5(I)-Y6(I))
C
         GRID_VZ(I,KM)=am(km)*rrho(km)*(Y11(I)-Y22(I))
         DIFF_NV(I,KM)=am(km)*rrho(km)*(Y33(I)-Y44(I))
         DIFF_VZ(I,KM)=am(km)*rrho(km)*(Y55(I)-Y66(I))
c
         y1(i)=y2(i)
         y3(i)=y4(i)
         y5(i)=y6(i)
         y11(i)=y22(i)
         y33(i)=y44(i)
         y55(i)=y66(i)
c
  200  continue
 2000 continue
c       a22k=am(kles)*rrho(kles)
      do 210 i=2,iles
       u1(i,kles)=am(kles)*rrho(kles)*(y2(i)+y4(i))
       v1(i,kles)=am(kles)*rrho(kles)*(y22(i)+y44(i))
C
c       GRID_UZ(I,KLES)=GRID_UZ(I,KLES)+am(kles)*rrho(kles)*Y2(i)
c       DIFF_NU(I,KLES)=DIFF_NU(I,KLES)+am(kles)*rrho(kles)*y4(i)
c       DIFF_UZ(I,KLES)=DIFF_UZ(I,KLES)+am(kles)*rrho(kles)*y6(i)
C
c       GRID_VZ(I,KLES)=GRID_VZ(I,KLES)+am(kles)*rrho(kles)*Y22(i)
c       DIFF_NV(I,KLES)=DIFF_NV(I,KLES)+am(kles)*rrho(kles)*y44(i)
c       DIFF_VZ(I,KLES)=DIFF_VZ(I,KLES)+am(kles)*rrho(kles)*y66(i)

       GRID_UZ(I,KLES)=am(kles)*rrho(kles)*Y2(i)
       DIFF_NU(I,KLES)=am(kles)*rrho(kles)*y4(i)
       DIFF_UZ(I,KLES)=am(kles)*rrho(kles)*y6(i)
C
       GRID_VZ(I,KLES)=am(kles)*rrho(kles)*Y22(i)
       DIFF_NV(I,KLES)=am(kles)*rrho(kles)*y44(i)
       DIFF_VZ(I,KLES)=am(kles)*rrho(kles)*y66(i)
C
  210 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-13-97)
        if (ijkadv .eq. 0 .and. I_FOUR .EQ. 1) then
c
c
               r2dt2=0.005/dt
c
           do i=2,iles
             do k=2,kles
               y1(k)=uu1(i,k+1)-2.*uu1(i,k)+uu1(i,k-1)
     1               -ub1(k+1)+2.*ub1(k)-ub1(k-1)
               y2(k)=vv1(i,k+1)-2.*vv1(i,k)+vv1(i,k-1)
     1               -vb1(k+1)+2.*vb1(k)-vb1(k-1)
             enddo
             y1(1)=y1(2)
             y2(1)=y2(2)
             y1(kmax)=y1(kles)
             y2(kmax)=y2(kles)
             do k=2,kles
                u1(i,k)=u1(i,k)-r2dt2*(y1(k+1)-2.*y1(k)+y1(k-1))
                v1(i,k)=v1(i,k)-r2dt2*(y2(k+1)-2.*y2(k)+y2(k-1))
c
                DIFF_NU(I,K)=DIFF_NU(I,K)
     1                         -R2DT2*(Y1(K+1)-2.*Y1(K)+Y1(K-1))
                DIFF_NV(I,K)=DIFF_NV(I,K)
     1                         -R2DT2*(Y2(K+1)-2.*Y2(K)+Y2(K-1))
c
             enddo
           enddo
cc        IF (I_FOUR .EQ. 1) THEN
c
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c
c           do i=2,iles
c             do k=3,kl2-1
c               y1(k)=(uu1(i,k+2)-3.*(uu1(i,k+1)-uu1(i,k))-uu1(i,k-1))
c     1               -(ub1(k+2)-3.*(ub1(k+1)-ub1(k))-ub1(k-1))
c               y2(k)=(vv1(i,k+2)-3.*(vv1(i,k+1)-vv1(i,k))-vv1(i,k-1))
c     1               -(vb1(k+2)-3.*(vb1(k+1)-vb1(k))-vb1(k-1))
c             enddo
c             do k=4,kl2-1
c                u1(i,k)=u1(i,k)+r2dt2*(y1(k-1)-y1(k))
c                v1(i,k)=v1(i,k)+r2dt2*(y2(k-1)-y2(k))
c             enddo
c             u1(i,3)=u1(i,3)+r2dt8*(uu1(i,4)-2.*uu1(i,3)+uu1(i,2))
c     1              -r2dt8*(ub1(4)-2.*ub1(3)+ub1(2))
c             v1(i,3)=v1(i,3)+r2dt8*(vv1(i,4)-2.*vv1(i,3)+vv1(i,2))
c             u1(i,kl2)=u1(i,kl2)+r2dt8*(uu1(i,kles)-2.*uu1(i,kl2)
c     1                          +uu1(i,kl2-1))
c     2                 -r2dt8*(ub1(kles)-2.*ub1(kl2)+ub1(kl2-1))
c
c             v1(i,kl2)=v1(i,kl2)+r2dt8*(vv1(i,kles)-2.*vv1(i,kl2)
c     1                          +vv1(i,kl2-1))
c     2                 -r2dt8*(vb1(kles)-2.*vb1(kl2)+vb1(kl2-1))
c          enddo

c
c        ELSE
c
c          y1(2)=0.
c          y1(kmax)=0.
c          y2(2)=0.
c          y2(kmax)=0.
c
c               r2dt8=0.008/d2t
c          do i=2,iles

c            do k=3,kles
c               y1(k)=uu1(i,k)-uu1(i,k-1)-ub1(k)+ub1(k-1)
c               y2(k)=vv1(i,k)-vv1(i,k-1)-vb1(k)+vb1(k-1)
c            enddo
c            do k=2,kles
c               u1(i,k)=u1(i,k)-r2dt8*(y1(k)-y1(k+1))
c               v1(i,k)=v1(i,k)-r2dt8*(y2(k)-y2(k+1))
c            enddo

c          enddo
C
C
c
cc        ENDIF
c
        endif
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nudge .ne. 1) then
        do 60 k=2,kles
          y1(k)=0.
   60     y3(k)=0.
        do 70 i=2,iles
        do 70 k=2,kles
          y1(k)=y1(k)+u1(i,k)
   70     y3(k)=y3(k)+v1(i,k)
        do 75 k=2,kles
          y1(k)=y1(k)*ril2
   75     y3(k)=y3(k)*ril2
        do 80 k=2,kles
        do 80 i=2,iles
          u1(i,k)=u1(i,k)-y1(k)
   80     v1(i,k)=v1(i,k)-y3(k)
      ENDIF
cc    ******   horizontal advections   *********************************
      do 2500 k=2,kles
       km=k-1
       kp=k+1
c
       y1(1)=a4k(k)*(uu1(3,k)-3.*(uu1(2,k)-uu1(1,k))-uu1(il2,k))
       y2(1)=a4k(k)*(vv1(3,k)-3.*(vv1(2,k)-vv1(1,k))-vv1(il2,k))
       do i=2,il2
         y1(i)=a4k(k)*(uu1(i+2,k)-3.*(uu1(i+1,k)-uu1(i,k))-uu1(i-1,k))
         y2(i)=a4k(k)*(vv1(i+2,k)-3.*(vv1(i+1,k)-vv1(i,k))-vv1(i-1,k))
       enddo
       do i=2,il2
         u1(i,k)=u1(i,k)+y1(i-1)-y1(i)
         v1(i,k)=v1(i,k)+y2(i-1)-y2(i)
c
          DIFF_NU(I,K)=DIFF_NU(I,K)+y1(i-1)-y1(i)
          DIFF_NV(I,K)=DIFF_NV(I,K)+y2(i-1)-y2(i)
c
       enddo
       u1(iles,k)=u1(iles,k)+y1(il2)-y1(1)
       v1(iles,k)=v1(iles,k)+y2(il2)-y2(1)
c
         DIFF_NU(ILES,K)=DIFF_NU(ILES,K)+Y1(IL2)-Y1(1)
         DIFF_NV(ILES,K)=DIFF_NV(ILES,K)+Y2(IL2)-Y2(1)
c
       do i=2,iles
          u1(i,k)=u1(i,k)-a_irf*rfa(k)*(uu1(i,k)-ub1(k))
     1                   +fcor*(v(i,k)-vb(k))
          v1(i,k)=v1(i,k)-a_irf*rfa(k)*(vv1(i,k)-vb1(k))
     2                   -fcor*(u(i,k)-ub(k))
c          u1(i,k)=u1(i,k)-a0000*am(k)*y0(i)
c     1          *(wb(kp)+wb(k))*(u(i,kp)-ub(kp)-u(i,km)+ub(km))*rd4z
c    2          -a0000*ub(k)*(u(i+1,k)-u(i-1,k))*rd2x
c          v1(i,k)=v1(i,k)-a0000*am(k)*y0(i)
c     1          *(wb(kp)+wb(k))*(v(i,kp)-vb(kp)-v(i,km)+vb(km))*rd4z
c     2          -a0000*ub(k)*(v(i+1,k)-v(i-1,k))*rd2x
c
          DIFF_NU(I,K)=DIFF_NU(I,K)-A_IRF*RFA(K)*(UU1(I,K)-UB1(K))
          DIFF_NV(I,K)=DIFF_NV(I,K)-A_IRF*RFA(K)*(VV1(I,K)-VB1(K))
c
       ENDDO
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       c33k=c3(k)*rdx
       ub2=0.0
c
       ub12=2.*ub(k)
       vb12=2.*vb(k)
c
       if (iadvh .eq. 2) then
       do 220 i=1,iles
c        y1(i)=(u(i,k)+u(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
c     1          -2.*ak(i,k)*(uu1(i+1,k)-uu1(i,k))*rdx2
c     2          +c33k*ak(i,k)*ak(i,k)
c        y2(i)=(v(i,k)+v(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
c     1          -(ak(i,k))*(vv1(i+1,k)-vv1(i,k))*rdx2
c
        y1(i)=(u(i,k)+u(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
        y3(i)=-2.*ak(i,k)*(uu1(i+1,k)-uu1(i,k))*rdx2
     1        +c33k*ak(i,k)*ak(i,k)
        y2(i)=(v(i,k)+v(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
        y4(i)=-(ak(i,k))*(vv1(i+1,k)-vv1(i,k))*rdx2
c
        y11(i)=(u(i,k)+u(i+1,k)-ub12)*(u(i,k)+u(i+1,k)-ub12)*rd4x
        y22(i)=(v(i,k)+v(i+1,k)-vb12)*(u(i,k)+u(i+1,k)-ub12)*rd4x
  220  continue
       do 230 i=2,iles
         u1(i,k)=u1(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i)
         v1(i,k)=v1(i,k)+y2(i-1)-y2(i)+y4(i-1)-y4(i)
c
c          GRID_UX(I,K)=GRID_UX(I,K)+Y1(I-1)-Y1(I)
c          DIFF_UX(I,K)=DIFF_UX(I,K)+Y11(I-1)-Y11(I)
c
c          GRID_VX(I,K)=GRID_VX(I,K)+Y2(I-1)-Y2(I)
c          DIFF_VX(I,K)=DIFF_VX(I,K)+Y22(I-1)-Y22(I)
c
c
          GRID_UX(I,K)=Y1(I-1)-Y1(I)
          DIFF_NU(I,K)=DIFF_NU(I,K)+Y3(I-1)-Y3(I)
          DIFF_UX(I,K)=Y11(I-1)-Y11(I)
c
          GRID_VX(I,K)=Y2(I-1)-Y2(I)
          DIFF_NV(I,K)=DIFF_NV(I,K)+Y4(I-1)-Y4(I)
          DIFF_VX(I,K)=Y22(I-1)-Y22(I)
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  230 continue
c
      else
c
cc    ******   4-th order horizontal advections   **********************
       do 240 i=1,iles
        y1(i)=-2.*ak(i,k)*(uu1(i+1,k)-uu1(i,k))*rdx2
     1         +c3(k)*ak(i,k)**2*rdx
        y2(i)=-ak(i,k)*(vv1(i+1,k)-vv1(i,k))*rdx2
  240  continue
       do 250 i=2,iles
        u1(i,k)=u1(i,k)+y1(i-1)-y1(i)
        v1(i,k)=v1(i,k)+y2(i-1)-y2(i)
c
          DIFF_NU(I,K)=DIFF_NU(I,K)+Y1(I-1)-Y1(I)
          DIFF_NV(I,K)=DIFF_NV(I,K)+Y2(I-1)-Y2(I)
c
  250  continue
       y1(1)=d516x*(u(1,k)+u(2,k))*(u(1,k)+u(2,k)-ub2)
     1       -d32x*(u(il2,k)+u(3,k))*(u(il2,k)+u(3,k)-ub2)
       y2(1)=d516x*(v(1,k)+v(2,k))*(u(1,k)+u(2,k)-ub2)
     1       -d32x*(v(il2,k)+v(3,k))*(u(il2,k)+u(3,k)-ub2)
c
       y11(1)=d516x*(u(1,k)+u(2,k)-ub12)*(u(1,k)+u(2,k)-ub12)
     1       -d32x*(u(il2,k)+u(3,k)-ub12)*(u(il2,k)+u(3,k)-ub12)
       y22(1)=d516x*(v(1,k)+v(2,k)-vb12)*(u(1,k)+u(2,k)-ub12)
     1       -d32x*(v(il2,k)+v(3,k)-vb12)*(u(il2,k)+u(3,k)-ub12)
c
       do 260 i=2,il2
        y1(i)=d516x*(u(i,k)+u(i+1,k))*(u(i,k)+u(i+1,k)-ub2)
     1        -d32x*(u(i-1,k)+u(i+2,k))*(u(i-1,k)+u(i+2,k)-ub2)
        y2(i)=d516x*(v(i,k)+v(i+1,k))*(u(i,k)+u(i+1,k)-ub2)
     1        -d32x*(v(i-1,k)+v(i+2,k))*(u(i-1,k)+u(i+2,k)-ub2)
c
        y11(i)=d516x*(u(i,k)+u(i+1,k)-ub12)*(u(i,k)+u(i+1,k)-ub12)
     1        -d32x*(u(i-1,k)+u(i+2,k)-ub12)*(u(i-1,k)+u(i+2,k)-ub12)
        y22(i)=d516x*(v(i,k)+v(i+1,k)-vb12)*(u(i,k)+u(i+1,k)-ub12)
     1        -d32x*(v(i-1,k)+v(i+2,k)-vb12)*(u(i-1,k)+u(i+2,k)-ub12)
  260  continue
       y3(1)=-d96x*(u(il2,k)+u(1,k))*(u(il2,k)+u(1,k)-ub2)
       y4(1)=-d96x*(v(il2,k)+v(1,k))*(u(il2,k)+u(1,k)-ub2)
c
       y33(1)=-d96x*(u(il2,k)+u(1,k)-ub12)*(u(il2,k)+u(1,k)-ub12)
       y44(1)=-d96x*(v(il2,k)+v(1,k)-vb12)*(u(il2,k)+u(1,k)-ub12)
       do 270 i=2,imax
        y3(i)=-d96x*(u(i-1,k)+u(i,k))*(u(i-1,k)+u(i,k)-ub2)
        y4(i)=-d96x*(v(i-1,k)+v(i,k))*(u(i-1,k)+u(i,k)-ub2)
c
        y33(i)=-d96x*(u(i-1,k)+u(i,k)-ub12)*(u(i-1,k)+u(i,k)-ub12)
        y44(i)=-d96x*(v(i-1,k)+v(i,k)-vb12)*(u(i-1,k)+u(i,k)-ub12)
  270  continue
       do 280 i=2,il2
        u1(i,k)=u1(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)
        v1(i,k)=v1(i,k)+y2(i-1)-y2(i)+y4(i-1)-y4(i+2)
c
c          GRID_UX(I,K)=GRID_UX(I,K)+Y1(I-1)-Y1(I)+Y3(I-1)-Y3(I+2)
c          GRID_VX(I,K)=GRID_VX(I,K)+Y2(I-1)-Y2(I)+Y4(I-1)-Y4(I+2)
c
c          DIFF_UX(I,K)=DIFF_UX(I,K)+Y11(I-1)-Y11(I)+Y33(I-1)-Y33(I+2)
c          DIFF_VX(I,K)=DIFF_VX(I,K)+Y22(I-1)-Y22(I)+Y44(I-1)-Y44(I+2)
c
          GRID_UX(I,K)=Y1(I-1)-Y1(I)+Y3(I-1)-Y3(I+2)
          GRID_VX(I,K)=Y2(I-1)-Y2(I)+Y4(I-1)-Y4(I+2)
c
          DIFF_UX(I,K)=Y11(I-1)-Y11(I)+Y33(I-1)-Y33(I+2)
          DIFF_VX(I,K)=Y22(I-1)-Y22(I)+Y44(I-1)-Y44(I+2)
  280  continue
       u1(iles,k)=u1(iles,k)+y1(il2)-y1(1)+y3(il2)-y3(3)
       v1(iles,k)=v1(iles,k)+y2(il2)-y2(1)+y4(il2)-y4(3)
c
c        GRID_UX(ILES,K)=GRID_UX(ILES,K)+Y1(IL2)-Y1(1)+Y3(IL2)-Y3(3)
c        GRID_VX(ILES,K)=GRID_VX(ILES,K)+Y2(IL2)-Y2(1)+Y4(IL2)-Y4(3)
c
c        DIFF_UX(ILES,K)=DIFF_UX(ILES,K)+Y11(IL2)-Y11(1)+Y33(IL2)-Y33(3)
c        DIFF_VX(ILES,K)=DIFF_VX(ILES,K)+Y22(IL2)-Y22(1)+Y44(IL2)-Y44(3)
c
        GRID_UX(ILES,K)=Y1(IL2)-Y1(1)+Y3(IL2)-Y3(3)
        GRID_VX(ILES,K)=Y2(IL2)-Y2(1)+Y4(IL2)-Y4(3)
c
        DIFF_UX(ILES,K)=Y11(IL2)-Y11(1)+Y33(IL2)-Y33(3)
        DIFF_VX(ILES,K)=Y22(IL2)-Y22(1)+Y44(IL2)-Y44(3)
       endif
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 2500 continue
c
c
ctao 9/22/97
       if (nudge .ne. 1) then
         do 295 k=2,kles
          do i=2,iles
c
c             U_LARGE(I,K)=U_LARGE(I,K)+UBT(K)
c             V_LARGE(I,K)=V_LARGE(I,K)+VBT(K)
             U_LARGE(I,K)=UBT(K)
             V_LARGE(I,K)=VBT(K)
c
             u1(i,k)=u1(i,k)+ubt(k)
             v1(i,k)=v1(i,k)+vbt(k)
          enddo
  295    continue
       endif
        if (nudge .eq. 1) then
          do k=2,kles
          do i=2,iles
c
c             U_LARGE(I,K)=U_LARGE(I,K)-FACTOR_NUDING*(UB(K)-UBI(K))
c             V_LARGE(I,K)=V_LARGE(I,K)-FACTOR_NUDING*(VB(K)-VBI(K))
c
            u1(i,k)=u1(i,k)-factor_nuding*(ub(k)-ubi(k))
            v1(i,k)=v1(i,k)-factor_nuding*(vb(k)-vbi(k))
c
             U_LARGE(I,K)=-FACTOR_NUDING*(UB(K)-UBI(K))
             V_LARGE(I,K)=-FACTOR_NUDING*(VB(K)-VBI(K))
         enddo
         enddo
       endif
ctao 9/22/97
c
c
cc    ******************************************************************
      do 300 k=2,kles
      do 300 i=2,iles
         VTEM_DT(i,k)=vv1(i,k)
       y1(i)=v(i,k)+eps*(-2.*v(i,k)+vv1(i,k))
       v(i,k)=vv1(i,k)+v1(i,k)*d2t
       vv1(i,k)=y1(i)+eps*v(i,k)
  300 continue
      do 400 k=1,kmax
      do 400 i=1,imax
  400   ak(i,k)=2.*ak(i,k)
      return
      end
