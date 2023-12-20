
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadv (x,u,w,itq)

CC    ****   COMPUTE ADVECTION OF DIFFERENT TRACERS
      integer nx,nz,itq
      PARAMETER (NX=514,NZ=43,ITT=244,NT=2880)

      real x(nx,nz),u(nx,nz),w(nx,nz)

      common/timestat/ ndt_stat,itime_ave,mmave

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      common/bxz/ imax,iles,il2,kmax,kles,kl2

      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg
      common/bcorr/ id

      real    qcl(nx,nz),qrn(nx,nz),qci(nx,nz),qcs(nx,nz),qcg(nx,nz)
      real    qcl1(nx,nz),qrn1(nx,nz),qci1(nx,nz),qcs1(nx,nz),
     1        qcg1(nx,nz),ww1(nx,nz)
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1i/ qci
      common/b1s/ qcs
      common/b1g/ qcg
      common/b2c/ qcl1
      common/b2r/ qrn1
      common/b2i/ qci1
      common/b2s/ qcs1
      common/b2g/ qcg1
      common/b2w/ ww1

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

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

      common/bls/ y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/bch/ it(nx),iv(nt),ics(nx,4),ibz(nx,4)

      real    snth(nz,4,7),sntv(nz,4,7),sntd(nz,4,7),snqh(nz,4,7),
     $ snqv(nz,4,7),snqd(nz,4,7),snhh(nz,4,7),snhv(nz,4,7),snhd(nz,4,7)
      common/bsts2/ snth,sntv,sntd,snqh,snqv,snqd,snhh,snhv,snhd

      common/bch1/ rby(7)

C
c      COMMON/Q_BUGT/ Q1_G_H(NX,NZ),Q1_G_V(NX,NZ),
c     1               Q1A_G_H(NX,NZ),Q1A_G_V(NX,NZ),
c     2               Q1_D_H(NX,NZ),Q1_D_V(NX,NZ),
c     3               Q1A_D_H(NX,NZ),Q1A_D_V(NX,NZ),
c     4               Q2_G_H(NX,NZ),Q2_G_V(NX,NZ),
c     5               Q2A_G_H(NX,NZ),Q2A_G_V(NX,NZ),
c     6               Q2_D_H(NX,NZ),Q2_D_V(NX,NZ),
c     7               Q2A_D_H(NX,NZ),Q2A_D_V(NX,NZ),
c     8               Q1_HYD(NX,NZ),Q2_HYD(NX,NZ),
c     9               Q1A_HYD(NX,NZ),Q2A_HYD(NX,NZ),
c     9               Q1_RAD(NX,NZ),Q1A_RAD(NX,NZ),
c     9               IBUDSEC,RBUD
C

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,kc,mt
      real    a0,a11,a22

      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC    ******   VERTICAL ADVECTION TERMS   ******************************

      a20=ndt_stat*ril2

      DO 1 I=2,ILES
    1  X4(I)=X(I,1)
      DO 1000 K=2,KLES
       A0=AM(K)*RDZ*RRHO(K)
      DO 10 I=2,IMAX
        X1(I)=U(I,K)
        X2(I)=W(I,K)
        X3(I)=W(I,K+1)
        X5(I)=X(I,K)
   10 CONTINUE
      DO 15 I=2,ILES
        Y1(I)=X(I+1,K)
        Y2(I)=X(I,K)
         IF (X1(I+1) .GE. 0.0) Y1(I)=X(I,K)
         IF (X1(I) .GT. 0.0) Y2(I)=X(I-1,K)
        Y3(I)=X(I,K+1)
        Y4(I)=X(I,K)
         IF (X3(I) .GE. 0.0) Y3(I)=X(I,K)
         IF (X2(I) .GT. 0.0) Y4(I)=X4(I)
        Y7(I)=-A0*(RHO1(K+1)*X3(I)*Y3(I)-RHO1(K)*X2(I)*Y4(I))
        Y8(I)=-RDX*(X1(I+1)*Y1(I)-X1(I)*Y2(I))
   15  CONTINUE

      DO I=2,ILES
c
c
c         if (ismg .eq. 1) then
c            q1_g_v(i,k)=q1_g_v(i,k)+y7(i)
c            q1_g_h(i,k)=q1_g_h(i,k)+y8(i)
c         endif
c         if (ismg .eq. 2) then
c            q2_g_v(i,k)=q2_g_v(i,k)+y7(i)
c            q2_g_h(i,k)=q2_g_h(i,k)+y8(i)
c         endif
c
c
         X(I,K)=X(I,K)+(Y7(I)+Y8(I))*DT
         X4(I)=X5(I)
         Y2(I)=.01*(W(I,K)+wb(k)*y0(i))
      ENDDO

      IF (ITQ .EQ. 1) THEN
       DO 20 I=2,ILES
         Y7(I)=-Y7(I)
         Y8(I)=-Y8(I)
   20  CONTINUE
      ENDIF
C     ******   STAT FOR W(QC+QR+QI+QS+QG)             ***************
c
       if(id.eq.1 .or. idq .eq. 1) then
        do 415 kc=1,7
         do 425 mt=1,4
         do 425 i=2,iles
            y16(i)=.01*(ww1(i,k)+wb(k)*y0(i))
          ibz(i,mt)=0
          ibz(i,1)=1
           if(ics(i,mt).eq.1) ibz(i,mt)=1
           if(kc.eq.4) go to 435
           if(kc.le.3) go to 445
            if (y16(i).gt.rby(kc)) ibz(i,mt)=0
            go to 435
  445       if (y16(i).lt.rby(kc)) ibz(i,mt)=0
  435     continue
  425     continue
c
         do 455 mt=1,4
          a11=0.
          a22=0.
         do 465 i=2,iles
          if(ibz(i,mt).eq.1) then
             a11=a11+y8(i)*a20
             a22=a22+y7(i)*a20
          endif
  465   continue
c
         if (ismg.eq.1) snth(k,mt,kc)=snth(k,mt,kc)+a11
         if (ismg.eq.1) sntv(k,mt,kc)=sntv(k,mt,kc)+a22
         if (ismg.eq.2) snqh(k,mt,kc)=snqh(k,mt,kc)+a11
         if (ismg.eq.2) snqv(k,mt,kc)=snqv(k,mt,kc)+a22
         if (ismg.gt.2) snhh(k,mt,kc)=snhh(k,mt,kc)+a11
         if (ismg.gt.2) snhv(k,mt,kc)=snhv(k,mt,kc)+a22

  455   continue
  415  continue
c
      endif
C     ******   end STAT FOR W(QC+QR+QI+QS+QG)             ***********
c
 1000 CONTINUE
      RETURN
      END
