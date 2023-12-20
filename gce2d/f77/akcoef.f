
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine akcoef (pi)
      parameter (NX=514,NZ=43,ITT=244)
      parameter (nb=nx*nz-6*nx)

      real    pi(nz)

      common/cpsbm/ icps_bm,iexplicit

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
    
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg

      common/bcor1/ fcor

      common/add_f/ i_four

      common/bstart/ n,isec,nran

      real    dpt1(nx,nz),dqv1(nx,nz),qcl1(nx,nz),qci1(nx,nz),
     $        uuT (nx,nz),wwT (nx,nz),ak1 (nx,nz),ak  (nx,nz)
      common/b2t/ dpt1
      common/b2q/ dqv1
      common/b2c/ qcl1
      common/b2i/ qci1
      common/b1a/ ak
      common/b2u/ uuT
      common/b2w/ wwT
      common/b2a/ ak1
   
      common/badv/ u1(nx,nz)

      common/badv1/ t0(nx),q0(nx),a0(nx),a(nx),a1(nx),a2(nx),ddd(nb)
   
      common/bsat/ uu1(nx,nz)
    
      common/bsat1/ ww1(nx,nz)
    
      common/dumuw/ umd(nx,nz),vmd(nx,nz),wmd(nx,nz)

      common/o4/ a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
   
      common/b4/ tbsk(nz),bskt(nz),bsk2t(nz),bsk(nz),bsk4(nz),
     1           bsit(nz),bsi2t(nz),bsi(nz),bsi4(nz)
   
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
    
      common/ba/ y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),t(nx),
     $        tp(nx),tm(nx),q(nx),qp(nx),qm(nx),y8(nx),y9(nx),y10(nx)
   
      common/bls/ y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)

      common/damp/ rfa(nz),rfa1(nz),cntc(nz),cgwd

c-----------------------------------------------------------------------

      save

c     ********************************
       a0000=0.
         if (IMLIFTING .eq. 1) a0000=1.
       a_irf=0.
         if (irf .eq. 1) a_irf=1.

      alsq=.622*al*al
      cpal=.622*1.61*cp*al
      alcp=al/cp
      cpra=cp*ra
      do 10 k=2,kles
       y1(k)=-am(k)*c1(k)*rd2z
   10  y3(k)=1./ta1(k)
      do 100 k=1,kmax
       y2(k)=1./pi(k)
      do 100 i=1,imax
  100  ak(i,k)=.5*ak(i,k)
c
       IF (IJKADV .EQ. 0) THEN
         do k=1,kmax
         do i=1,imax
            UU1(I,K)=UUT(I,K)
            WW1(I,K)=WWT(I,K)
         enddo
         enddo
       ELSE
         do k=1,kmax
         do i=1,imax
            UU1(I,K)=.5*(3.*UUT(I,K)-UMD(I,K))
            WW1(I,K)=.5*(3.*WWT(I,K)-WMD(I,K))
         enddo
         enddo
       ENDIF
c
c     ********************************
      do 210 i=2,iles
       tm(i)=ta1(1)+dpt1(i,1)
       t(i)=ta1(2)+dpt1(i,2)
       qm(i)=qa1(1)+dqv1(i,1)
       q(i)=qa1(2)+dqv1(i,2)
  210 continue
      do 2000 k=2,kles
        kp=k+1
        km=k-1
c        y13k=y1(k)*y3(k)
       do 250 i=2,iles
         y8(i)=qcl1(i,k)+qci1(i,k)
        tp(i)=ta1(kp)+dpt1(i,kp)
        qp(i)=qa1(kp)+dqv1(i,kp)
        if (y8(i).ge.1.e-5) then
          t0(i)=t(i)*pi(k)
          q0(i)=q(i)
          a0(i)=cpra*t0(i)*t0(i)
          a(i)=(a0(i)+cpal*t0(i)*q0(i))/((a0(i)+alsq*q0(i))*ta1(k))
          a1(i)=tp(i)-tm(i)+alcp*(qp(i)*y2(kp)-qm(i)*y2(km))
          a2(i)=qp(i)-qm(i)+qcl1(i,kp)-qcl1(i,km)+qci1(i,kp)-qci1(i,km)
         u1(i,k)=y1(k)*(a(i)*a1(i)-a2(i))
        else
         u1(i,k)=y1(k)*(tp(i)*(1.+.61*qp(i))-tm(i)*(1.+.61*qm(i)))*y3(k)
        endif
         tm(i)=t(i)
         t(i)=tp(i)
         qm(i)=q(i)
         q(i)=qp(i)
  250  continue
 2000 continue

      do 310 i=1,imax
       qm(i)=ak1(i,1)*ak1(i,1)
       tm(i)=ak1(i,2)*ak1(i,2)
  310 continue

      do 3000 k=2,kles
        kp=k+1
        km=k-1
        a11k=am(k)*rd4z
c        c11k=2.*coef(k)*am(k)*am(k)*rdz2
       do 320 i=1,imax
  320   tp(i)=ak1(i,kp)*ak1(i,kp)
       do 330 i=2,iles
         y9(i)=a0000*y0(i)*wb(k)
         y10(i)=a0000*y0(i)*wb(kp)
       y1(i)=-c2(k)*tm(i)+coef(k)*
     1   (((ww1(i+1,kp)+ww1(i+1,k)-ww1(i-1,kp)-ww1(i-1,k))*rd4x
     2     +(uu1(i+1,kp)+uu1(i,kp)-uu1(i+1,km)-uu1(i,km))*a11k)**2
     3    +2.*(uu1(i+1,k)-uu1(i,k))**2*rdx2)
     4   +(tm(i+1)+tm(i-1)-2.*tm(i))*r2dx2
     5   +am(k)*(am1(kp)*(tp(i)-tm(i))-am1(k)*(tm(i)-qm(i)))*r2dz2
     6   +2.*coef(k)*am(k)*am(k)*rdz2
     7                 *(ww1(i,kp)-ww1(i,k)+y10(i)-y9(i))**2
       u1(i,k)=u1(i,k)+y1(i)
       u1(i,k)=u1(i,k)-a0000*(a11k*y0(i)*
     1                   (wb(k+1)+wb(k))*(ak(i,k+1)-ak(i,k-1)))
  330  continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       do 390 i=1,imax
        qm(i)=tm(i)
  390   tm(i)=tp(i)
 3000 continue
c
      if (ijkadv .eq. 1) then
c
       awks=.5
c
      do 5000 k=2,kles
        kp=k+1
        km=k-1
c         y9(i)=y0(i)*wb(k)
c         y10(i)=y0(i)*wb(kp)
c          IF (IMLIFTING .EQ. 0) THEN
           do i=2,iles
             y9(i)=0.0
             y10(i)=0.0
           enddo
c          ENDIF
        a22z=AM(K)*RRHO(K)*RD2Z
        ammkp=awks*am1(kp)*rho1(kp)*bsk(kp)
        ammk=awks*am1(k)*rho1(k)*bsk(k)
       do 520 i=2,iles
        u1(i,k)=u1(i,k)-a22z*(
     1   rho1(kp)*y10(i)*(ak(i,k)+ak(i,kp))
     2   -rho1(k)*y9(i)*(ak(i,k)+ak(i,km))
     3   -2.*(ammkp*(ak1(i,kp)-ak1(i,k))-ammk*(ak1(i,k)-ak1(i,km)))*rdz)
        U1(I,K)=U1(I,K)-a_irf*RFA(K)*AK1(I,K)
  520  continue
        y2(1)=a4k(k)*(ak1(3,k)-3.*(ak1(2,k)-ak1(1,k))-ak1(il2,k))
        do 540 i=2,il2
         y2(i)=a4k(k)*(ak1(i+2,k)-3.*(ak1(i+1,k)-ak1(i,k))-ak1(i-1,k))
  540   continue
        do 570 i=2,il2
  570    u1(i,k)=u1(i,k)+awks*(y2(i-1)-y2(i))
        u1(iles,k)=u1(iles,k)+awks*(y2(il2)-y2(1))

        do i=2,imax
           y1(i)=-(ak(i,k)+ak(i-1,k))*(ak1(i,k)-ak1(i-1,k))*r2dx2
        enddo
        do i=2,iles
           u1(i,k)=u1(i,k)+y1(i)-y1(i+1)
        enddo

 5000 continue
c
       do i=2,iles
          y1(i)=0.
       enddo
      do 4000 k=3,kmax
        km=k-1
       do i=2,iles
          y2(i)=rho1(k)*(-(am1(k)*(ak1(i,k)-ak1(i,km))
     1                  *(ak(i,k)+ak(i,km)))*rdz)
       enddo
       if (km .eq. kles) then
          do i=2,iles
             y2(i)=0.
          enddo
       endif
       do i=2,iles
        u1(i,km)=u1(i,km)+am(km)*(y1(i)-y2(i))*rd2z*rrho(km)
         y1(i)=y2(i)
       enddo
 4000 continue
cc
       call advectak (ak,ak1)
cc
      else
cc
      DO 6000 K=2,KLES
        KP=K+1
        KM=K-1
c        a22z=AM(K)*RRHO(K)*RD2Z
        ammkp=AM1(KP)*RHO1(KP)*BSK(KP)
        ammk=AM1(K)*RHO1(K)*BSK(K)
c         y9(i)=y0(i)*wb(k)
c         y10(i)=y0(i)*wb(kp)
c        IF (IMLIFTING .EQ. 0) THEN
           do i=2,iles
             y9(i)=0.0
             y10(i)=0.0
           enddo
c        ENDIF
       DO 610 I=2,ILES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          u1(i,k)=u1(i,k)-am(k)*rrho(k)*rd2z*(
     1      rho1(kp)*(ww1(i,kp)+y10(i))*(ak(i,k)+ak(i,kp))-2.
     2      *(ammkp*(ak1(i,kp)-ak1(i,k))-ammk*(ak1(i,k)-ak1(i,km)))*rdz        
     3      -rho1(k)*(ww1(i,k)+y9(i))*(ak(i,k)+ak(i,km)))
          U1(I,K)=U1(I,K)-a_irf*RFA(K)*AK1(I,K)
  610  CONTINUE
       if(iadvh .eq. 2) then
        do 620 i=2,iles
         u1(i,k)=u1(i,k)
     1   -(uu1(i+1,k)*(ak(i+1,k)+ak(i,k))-uu1(i,k)*(ak(i,k)+ak(i-1,k)))
c    2   *rd2x+bsi(k)*rdx2*(ak1(i+1,k)-2.*ak1(i,k)+ak1(i-1,k))
     2   *rd2x
  620   continue
        y2(1)=a4k(k)*(ak1(3,k)-3.*(ak1(2,k)-ak1(1,k))-ak1(il2,k))
        do 615 i=2,il2
         y2(i)=a4k(k)*(ak1(i+2,k)-3.*(ak1(i+1,k)-ak1(i,k))-ak1(i-1,k))
  615   continue
        do 625 i=2,il2
  625    u1(i,k)=u1(i,k)+y2(i-1)-y2(i)
        u1(iles,k)=u1(iles,k)+y2(il2)-y2(1)
       else
c     ***   4-th order horizontal advection scheme   ************
        y2(1)=uu1(2,k)*(d58x*(ak(2,k)+ak(1,k))-d16x*(ak(3,k)+ak(il2,k)))
     1        +a4k(k)*(ak1(3,k)-3.*(ak1(2,k)-ak1(1,k))-ak1(il2,k))
        do 630 i=2,il2
         y2(i)=uu1(i+1,k)*(d58x*(ak(i+1,k)+ak(i,k))
     1                     -d16x*(ak(i+2,k)+ak(i-1,k)))
     2         +a4k(k)*(ak1(i+2,k)-3.*(ak1(i+1,k)-ak1(i,k))-ak1(i-1,k))
  630   continue
        y3(1)=-d48x*uu1(1,k)*(ak(1,k)+ak(il2,k))
        do 640 i=2,imax
  640    y3(i)=-d48x*uu1(i,k)*(ak(i,k)+ak(i-1,k))
        do 650 i=2,il2
  650    u1(i,k)=u1(i,k)+y2(i-1)-y2(i)+y3(i-1)-y3(i+2)
        u1(iles,k)=u1(iles,k)+y2(il2)-y2(1)+y3(il2)-y3(3)
       endif
 6000 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-13-97)

c        IF (I_FOUR .EQ. 1) THEN
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c           do i=2,iles
c             do k=3,kl2-1
c               y1(k)=(ak1(i,k+2)-3.*(ak1(i,k+1)-ak1(i,k))-ak1(i,k-1))
c             enddo
c             do k=4,kl2-1
c                u1(i,k)=u1(i,k)+r2dt2*(y1(k-1)-y1(k))/d2t
c             enddo
c             u1(i,3)=u1(i,3)+r2dt8*(ak1(i,4)-2.*ak1(i,3)+ak1(i,2))
c             u1(i,kl2)=u1(i,kl2)+r2dt8*(ak1(i,kles)-2.*ak1(i,kl2)
c     1                                 +ak1(i,kl2-1))
c          enddo
c
c
c        ELSE
c 
c          y1(2)=0.
c          y1(kmax)=0.

c               r2dt8=0.008/d2t

c            do i=2,iles

c            do k=3,kles
c               y1(k)=ak1(i,k)-ak1(i,k-1)
c            enddo
c            do k=2,kles
c               u1(i,k)=u1(i,k)-r2dt8*(y1(k)-y1(k+1))
c            enddo

c          enddo
C
C
cc        ENDIF
c
        IF (I_FOUR .EQ. 1) THEN
           R2DT2=0.005/DT
           DO I=2,ILES
             DO K=2,KLES
               Y1(K)=AK1(I,K+1)-2.*AK1(I,K)+AK1(I,K-1)
             ENDDO
             Y1(1)=Y1(2)
             Y1(KMAX)=Y1(KLES)
             DO K=2,KLES
                U1(I,K)=U1(I,K)-R2DT2*(Y1(K+1)-2.*Y1(K)+Y1(K-1))
             ENDDO
           ENDDO
        ENDIF
cc
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      endif
c
c     ***********************************************************
      do 400 k=2,kles
       a11=tbsk(k)
      do 400 i=2,iles
       IF (IJKADV .EQ. 0) THEN
         Y1(I)=AK(I,K)+EPS*(-2.*AK(I,K)+AK1(I,K))
        AK(I,K)=AK1(I,K)+U1(I,K)*D2T
        AK1(I,K)=Y1(I)+EPS*AK(I,K)
       ELSE
        AK(I,K)=AK(I,K)+U1(I,K)*DT
c        ak1(i,k)=ak(i,k)
       ENDIF
        ak(i,k)=min(a11, max(ak(i,k), 0.0))
        ak1(i,k)=min(a11, max(ak1(i,k), 0.0))
  400 continue
      call boundy (ak,ak1)
      return
      end
