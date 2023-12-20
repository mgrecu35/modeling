cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine opt4(i1,nx1,nx2,pl,ta,wa,oa,tausw,tauir
     +               ,fcld,taq,waq,oaq,plq,tsq,reff)
c define variables and calculate the optical thickness
      integer nx,nz,lay,nadd,nz2,nz4,nz15,nx3
      integer npp1,mb1,mb2,nq,nw
      parameter(NX=514,NZ=43,lay=88)
      parameter(nadd=7,nz2=nz*2,nz4=nz*4,nz15=nz*15,nx3=nx*3)
      parameter(npp1=nz+nadd-1,mb1=nz15,mb2=nz*7+nx3,nq=npp1,nw=nq-1)

      integer i1,nx1,nx2
      real    pl(lay),ta(lay),wa(lay),oa(lay)
      real    tausw(nx1,1,nw,2),tauir(nx1,1,nw),reff(nx1,1,nw,2)
      real    fcld(nx1,1,nw),tsq(nx1,1)
      real    plq(nx1,1,nq),taq(nx1,1,nw),waq(nx1,1,nw),oaq(nx1,1,nw)

      integer iradave
      common/iptionr/ iradave
      COMMON/IPTIONR1/ IOPCLOUD
      common/iceopt/ ice913,ilif
      real    tnw,tns,tng,roqs,roqg,roqr
      common/size/ tnw,tns,tng,roqs,roqg,roqr

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),
     1        thairsf(nx),pairsfc(nx)
      common/dinrad/ p00,dz0,tairsfc,qairsfc,thairsf,pairsfc
      real    dpt(nx,nz),dqv(nx,nz),qci(nx,nz),qcl(nx,nz),
     $        qrn(nx,nz),qcs(nx,nz),qcg(nx,nz)
      common/b1t/ dpt
      common/b1q/ dqv
      common/b1i/ qci
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1s/ qcs
      common/b1g/ qcg
      real    tb(nz),qb(nz),zz0(nz),rho(nz),zz1(nz2),tb1(nz),
     $        qb1(nz),zz2(mb1)
      common/b5/ tb,qb,zz0,rho,zz1,tb1,qb1,zz2
      real    zz3(nz4),p0(nz),pi(nz),zz4(mb2)
      common/b6/ zz3,p0,pi,zz4
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
c     real    qci(nx1,nz),qcl(nx1,nz),qrn(nx1,nz),qcs(nx1,nz)
c     real    aa1(nx,npp1),bb(nx,npp1),qcg(nx1,nz)
      real    aa1(nx,npp1),bb(nx,npp1),aa2(nx,npp1)
      real    y1(nw),y2(nw),y3(nw),y4(nw),y5(nw),y6(nw),y7(nw)
      integer i,i2,i3,ii,k,km
      real    b0,b1,b2,b4,b5,b6,b7,b9,cpi,effrad
      real    q1,q2,q3,q5,rnx1,tauqc,tauqg,tauqii,tauqis,tauqr,twco
      save
c
c      twcz=1.e-5
      I2=I1+NX2-1
      DO I=I1,I2,2
        II=(I-I1)/2+1
        I3=I+1
        do k=1,npp1
          plq(ii,1,k)=pl(k)
        enddo
        do k=1,nw
          oaq(ii,1,k)=oa(k)
          fcld(ii,1,k)=1.0
        enddo
        do k=1,nadd
          waq(ii,1,k)=wa(k)
          taq(ii,1,k)=ta(k)
c          taq(ii,1,k)=190.0
        enddo
        do k=2,kles
           km=kmax-k+nadd
          WAQ(II,1,KM)=DQV(I3,K)+QB1(K)
            IF(WAQ(II,1,KM).LE.1.E-6) WAQ(II,1,KM)=1.E-6
          TAQ(II,1,KM)=(DPT(I3,K)+TB1(K))*PI(K)
            IF(TAQ(II,1,KM).LE.190.) TAQ(II,1,KM)=190.
        enddo
        tsq(ii,1)=tairsfc(i3)
         if (ilif .eq. 1) tsq(ii,1)=29.+273.16
      enddo
       cpi =4.*atan(1.)
c       twcz=1.e-5
       twco=1.e-6
       DO 200 I=1,NX1
          II=I1+(I-1)*2
          i3=ii+1
         do k=1,nw
           tausw(i,1,k,1)=0.
           tausw(i,1,k,2)=0.
           reff(i,1,k,1)=0.
           reff(i,1,k,2)=0.
           tauir(i,1,k)=0.
         enddo
c        do k=2,kles
c        if (tqe(i,k) .ge. twcz) go to 30
c        enddo
c        go to 200
c30      continue
         b5=0.
         b6=0.
         b7=0.
         b9=0.
c         ttauqc=0.
c         ttauqr=0.
c         ttauis=0.
c         ttauii=0.
c         ttauqg=0.
         do 100 k=2,kles
           km=kmax-k+nadd
c         tausw(i,1,km)=0.
c         tauir(i,1,km)=0.
            b0=0.
            b1=0.
            b2=0.
            b4=0.
            b3=0.
            tauqc=0.0
            tauqr=0.0
            tauqs=0.0
            tauqis=0.0
            tauqii=0.0
            tauqg=0.0
           q1=qcl(i3,k)
           q2=qrn(i3,k)
           IF (IOPCLOUD .EQ. 1) THEN
              Q3=QCI(I3,K)+QCS(I3,K)
              Q4=0.
           ELSE
              Q3=QCI(I3,K)
              Q4=QCS(I3,K)
           ENDIF
           q5=qcg(i3,k)
          if(q1 .ge. twco) then
            b0=rho(k)*dz0(k)*q1
           effrad=0.0015
           tauqc=b0/effrad
            b0=1.e4*b0
            reff(i,1,km,2)=effrad*1.0e4
          endif
          if(q2 .ge. twco) then
            b1=rho(k)*dz0(k)*q2
           effrad=3./((cpi*tnw*roqr/(rho(k)*q2))**.25)
           tauqr=b1/effrad
            b1=1.e4*b1
          endif
          if(q4 .ge. twco) then
            b3=rho(k)*dz0(k)*q4
           effrad=3./((cpi*tns*roqs/(rho(k)*q4))**.25)
           tauqs=b3/effrad
            b3=1.e4*b3
          endif
          if(q5 .ge. twco) then
            b4=rho(k)*dz0(k)*q5
           effrad=3./((cpi*tng*roqg/(rho(k)*q5))**.25)
           tauqg=b4/effrad
            b4=1.e4*b4
          endif
            b4=b4+b3
            tauqg=tauqg+tauqs
          if(q3 .ge. twco) then
            b2=1.e4*rho(k)*dz0(k)*q3
c          effrad=0.0050
           effrad=0.0125+(taq(i,1,km)-243.16)*0.00050
           if (taq(i,1,km) .gt. 243.16) effrad=0.0125
           if (taq(i,1,km) .lt. 223.16) effrad=0.0025
           tauqis=b2*(-0.006656+ 3.686e-4/effrad)
           tauqii=b2*(-0.011500+ 4.110e-4/effrad
     +                         +17.300e-8/(effrad*effrad))
          reff(i,1,km,1)=effrad*1.0e4
          endif
          b5=b5+b0
          b6=b6+b1
          b7=b7+b2
          b9=b9+b4
c           ttauqc=ttauqc+tauqc
c           ttauqr=ttauqr+tauqr
c           ttauis=ttauis+tauqis
c           ttauii=ttauii+tauqii
c           ttauqg=ttauqg+tauqg
c          tausw(i,1,km)=1.5*(tauqc+tauqr)
           tausw(i,1,km,2)=1.5*(tauqc+tauqr)
           tauir(i,1,km)=0.5*tausw(i,1,km,2)
c          tauir(i,1,km)=0.5*tausw(i,1,km)
c          tausw(i,1,km)=tausw(i,1,km)+tauqis+tauqg
           tausw(i,1,km,1)=tauqis+tauqg
           tauir(i,1,km)=tauir(i,1,km)+tauqii+tauqg
  100    continue
  200    continue
         DO I=1,NX1
           k=2
           km=kmax-k+nadd
           k=3
           km1=kmax-k+nadd
           aa1(i,km)=0.8*tausw(i,1,km,1)+0.2*tausw(i,1,km1,1)
           aa2(i,km)=0.8*tausw(i,1,km,2)+0.2*tausw(i,1,km1,2)
           bb(i,km)=0.8*tauir(i,1,km)+0.2*tauir(i,1,km1)
cc           aa1(i,km)=tausw(i,1,km,1)
cc           aa2(i,km)=tausw(i,1,km,2)
cc           bb(i,km)=tauir(i,1,km)
           do 110 k=3,kles-1
             km=kmax-k+nadd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            aa1(i,km)=0.2*tausw(i,1,km-1,1)+0.60*tausw(i,1,km,1)
     1                +0.2*tausw(i,1,km+1,1)
             aa2(i,km)=0.2*tausw(i,1,km-1,2)+0.60*tausw(i,1,km,2)
     1                +0.2*tausw(i,1,km+1,2)
             bb(i,km)=0.2*tauir(i,1,km-1)+0.60*tauir(i,1,km)
     1               +0.2*tauir(i,1,km+1)
cc             aa1(i,km)=tausw(i,1,km,1)
cc             aa2(i,km)=tausw(i,1,km,2)
cc             bb(i,km)=tauir(i,1,km)
  110      continue
           k=kles
           km=kmax-k+nadd
           k=kl2
           km1=kmax-k+nadd
            aa1(i,km)=0.8*tausw(i,1,km,1)+0.2*tausw(i,1,km1,1)
            aa2(i,km)=0.8*tausw(i,1,km,2)+0.2*tausw(i,1,km1,2)
            bb(i,km)=0.8*tauir(i,1,km)+0.2*tauir(i,1,km1)
cc           aa1(i,km)=tausw(i,1,km,1)
cc           aa2(i,km)=tausw(i,1,km,2)
cc           bb(i,km)=tauir(i,1,km)
           do 120 k=2,kles
              km=kmax-k+nadd
             tausw(i,1,km,1)=aa1(i,km)
             tausw(i,1,km,2)=aa2(i,km)
             tauir(i,1,km)=bb(i,km)
  120      continue

c          do k=1,nadd
c            tausw(i,1,k)=0.0
c            tauir(i,1,k)=0.0
c          enddo
        enddo
      IF (IRADAVE .EQ. 1) THEN
         RNX1=1./FLOAT(NX1)
        DO K=1,NW
          Y1(K)=0.
          Y2(K)=0.
          Y3(K)=0.
          Y4(K)=0.
          Y5(K)=0.
          Y6(K)=0.
          Y7(K)=0.
        ENDDO
        DO 300 K=2,KLES
          KM=KMAX-K+NADD
          DO I=1,NX1
             Y1(KM)=Y1(KM)+TAQ(I,1,KM)
             Y2(KM)=Y2(KM)+WAQ(I,1,KM)
             Y3(KM)=Y3(KM)+TAUSW(I,1,KM,1)
             Y4(KM)=Y4(KM)+TAUSW(I,1,KM,2)
             Y5(KM)=Y5(KM)+REFF(I,1,KM,1)
             Y6(KM)=Y6(KM)+REFF(I,1,KM,2)
             Y7(KM)=Y7(KM)+TAUIR(I,1,KM)
          ENDDO
  300   CONTINUE
         DO K=2,KLES
          KM=KMAX-K+NADD
          TAQ(1,1,KM)=Y1(KM)*RNX1
          WAQ(1,1,KM)=Y2(KM)*RNX1
          TAUSW(1,1,KM,1)=Y3(KM)*RNX1
          TAUSW(1,1,KM,2)=Y4(KM)*RNX1
          REFF(1,1,KM,1)=Y5(KM)*RNX1
          REFF(1,1,KM,2)=Y6(KM)*RNX1
          TAUIR(1,1,KM)=Y7(KM)*RNX1
        ENDDO
      ENDIF
      return
      end
