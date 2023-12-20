C-----------------------------------------------------------------------------
      subroutine base (IRS)
      PARAMETER (NX=514,NZ=43,ITT=244,NX10=10*NX)
      common/msui/ isounding,isui,ifine,idelz
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      common/iceopt/ ice913,ilif
      common/cpsbm/ icps_bm,iexplicit
      common/o4/ a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/b1u/ u(nx,nz)
      common/b1v/ v(nx,nz)
      common/b1w/ w(nx,nz)
      common/b2u/ uu1(nx,nz)
      common/b2v/ vv1(nx,nz)
      common/b2w/ ww1(nx,nz)
      common/dumuw/ umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/b4/ tbsk(nz),bskt(nz),bskt2(nz),bskm(nz),bskm4(nz),
     $   bsit(nz),bsit2(nz),bsim(nz),bsim4(nz)
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     $   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/bb6/ tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     $  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     $  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     *  q2t(nz),vb_6h(nz,itt),vbt(nz)

      COMMON/UPPER1/ T_ADJUST(NZ,ITT),Q_ADJUST(NZ,ITT),PRESS_T(NZ),
     1     TEMP_T(NZ),VAP_T(NZ)
      common/bls3/ factor_nuding
      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/q1q2z/ Q1Z_6H(NZ,ITT),Q1ZT(NZ),Q2Z_6H(NZ,ITT),Q2ZT(NZ)
      common/ba/ y1(nx),y2(nx),y3(nx),qvs(nx),qr(nx),hu1(nx),y4(nx10)
      common/damp/ rfa(nz),rfa1(nz),cntc(nz),cgwd
      common/bls/ wbx(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/dinrad/p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),thairsf(nx)
     *,       pairsfc(nx)
      common/dinrad1/ dz1half
      common/bpbl/ UHT(NZ),WHT(NZ),TGBAT0
      real z1(nz),z2(nz),hgt1(nz),hgt2(nz)
      common/zlevel/ cc1,cc2,z1,z2
      common/bnumer/ kref
      real dzy(nz),dzyp(nz),sst(itt),ssp(itt)
      save
c      DATA DZY/0.,120.,150.,180.,210.,240.,270.,300.,350.,400.,450.,
c     1   500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,
c     2   1000.,1000.,1000.,1000.,1000.,1100.,1100.,1100.,1100.,1100.,
c     3   1200.,1200./
c      DATA DZYP/0.,120.,150.,180.,210.,240.,270.,300.,350.,400.,450.,
c     1   500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,
c     2   1000.,1000.,1000.,1000.,1000.,1100.,1100.,1100.,1100.,1100.,
c     3   1200.,1200./

      ILES=NX-2
      KLES=NZ-2
      IL2=ILES
      KL2=KLES
      ILES=ILES+1
      KLES=KLES+1
      IMAX=ILES+1
      KMAX=KLES+1
      RIL2=1./FLOAT(IL2)

C Grid setup for GCSS WG4 model intercomparison workshop

      dzy(1)=0.
      dzyp(1)=0.
      dzy(2)=0.
      dzyp(2)=0.
      do k=3,nz
         if(k.le.22)then
            dzk=real(k-2)*300.
            dzkp=real(k-2.5)*300.
            dzy(k)=(.2982456142+1.169590643E-4*dzk)*dzk
            dzyp(k)=(.2982456142+1.169590643E-4*dzkp)*dzkp
         else
            dzy(k)=real(k-22)*500.+6000.
            dzyp(k)=real(k-22.5)*500.+6000.
         endif
      enddo
      do k=2,nz-1
         dzyp(k)=dzyp(k+1)
      enddo
      dzyp(nz)=dzyp(nz-1)+500.

C End GCSS WG4 grid setup

      if(idelz.eq.1)then 
         DZ=DZY(3)*100.
      else

C  IF USING 34 PTS, SET DZ = 27.e2
C  IF USING 40 PTS, SET DZ = 25.e2
C  IF USING 43 PTS, SET DZ = 20.e2

         DZ=20.e2      !cm

      endif
cccshie
      DX=1000.e2    !cm for nx= 514, 5/18-5/26/98
c     DX= 500.e2    !cm for nx= 1026, 6/2-6/11/98
c     DX=2000.e2    !cm for nx= 258, 5/18-5/20/98
c
      IF(ICPS_BM .EQ. 1) DX=30000.E2
c
      PSFC=1007.e3
      TSFC=27.2     !C

C IF USING 34 PTS, SET CC1=4, CC2=1/40.
C IF USING 40 PTS, SET CC1=6, CC2=1/50.
C IF USING 43 PTS, SET CC1=3.4, CC2=.03

      CC1=3.4
      CC2=.03

      RD1=1.E-3
      RD2=2.2
      BOUND=1.E-3
      AL=2.5E10
      CP=1.004E7
      RA=2.87E6
      F5=17.26939*237.3*AL/CP
      SCALE=-1.*1000./3600.

      Z1(1) = 0.
      Z2(1) = 0.

      if (idelz.eq.1) then
         do k=1,kmax
            Z2(k)=dzy(k)*100.
            Z1(k)=dzyp(k)*100.
         enddo
      else
         DO K=2,KMAX
            if (ifine .eq. 0) then
               Z1(k)=(cc1+cc2*(k-1.5)*dz*.01)*(k-1.5)*dz
               Z2(k)=(cc1+cc2*(k-2.)*dz*.01)*(k-2.)*dz
            else
               y1(k)=(cc1+cc2*(k-1.5)*dz*.01)*(k-1.5)*dz
               y2(k)=(cc1+cc2*(k-2.)*dz*.01)*(k-2.)*dz
            endif

         ENDDO
         if (ifine .eq. 1) then
            Z1(2)=40.E2
            Z2(2)=0.0
            Z2(3)=80.E2
            DO K=2,KLES
               Z1(K+1)=Y1(K)+80.E2
               Z2(K+1)=Y2(K)+80.E2
               print*,k,z2(k)/100.,(z2(k)-z2(k-1))/100.
            ENDDO
         endif
      endif

      DO K=2,KLES
         AM(K)=DZ/(Z2(K+1)-Z2(K))
         AM1(K)=DZ/(Z1(K)-Z1(K-1))
      ENDDO
      AM1(KMAX)=DZ/(Z1(KMAX)-Z1(KLES))
      AM1(1)=AM1(2)
      AM(1)=AM(2)
      AM(KMAX)=AM(KLES)
      do k=1,kmax
        hgt1(k)=z1(k)
        hgt2(k)=z2(k)
        uht(k)=z1(k)
        wht(k)=z2(k)
      enddo
      dz1half=z1(2)/100.
      write(*,*) z2
      write(6,*)' km      am        am1        z1          z2       
     1dz2       dratio'
      do 14 k=kmax,2,-1
      deltz2=z2(k)-z2(k-1)
      if(k.ge.4)then
         deltz3=z2(k-1)-z2(k-2)
         dratio=deltz2/deltz3
      else
         dratio=0.
      endif
   14 write(6,2035) k,am(k),am1(k),z1(k)/100.,z2(k)/100.,
     1              deltz2/100.,dratio
      write(6,*)
      write(6,*)

      DX2=DX*DX
      DZ2=DZ*DZ
      RDX=1./DX
      RDZ=1./DZ
      RD2X=.5*RDX
      RD2Z=.5*RDZ
      RD4X=.25*RDX
      RD4Z=.25*RDZ
      RDX2=RDX*RDX
      RDZ2=RDZ*RDZ
      R2DX2=.5*RDX2
      R2DZ2=.5*RDZ2
      R4DX2=.25*RDX2
      R4DZ2=.25*RDZ2
      D16X=1./16.*RDX
      D48X=1./48.*RDX
      D58X=5./8.*RDX
      D32X=1./32.*RDX
      D96X=1./96.*RDX
      D516X=5./16.*RDX

C     ******   RAYLEIGH RELAXATION TERMS 
      RRTP0=10.
      RRTP=15.5  ! tao 5/1/01
c     RRTP=16.
      RRTPD=RRTP-RRTP0
      CGWD=1./(5.*24.*3600.)
      CG=1./(24.*3600.)
      CRF=1./900.
      DO 15 K=1,KMAX
        A1=1.E-5*Z1(K)
        A2=1.E-5*Z2(K)
        UHT(K)=.01*Z1(K)
        WHT(K)=.01*Z2(K)
        RFA(K)=0.0
        RFA1(K)=0.0
        CNTC(K)=0.0
      IF(A1.GE.RRTP) CNTC(K)=CG
      IF(A1.LT.RRTP .AND. A1.GT.RRTP0) CNTC(K)=CG*(A1-RRTP0)/RRTPD
      IF(A1.GE.RRTP) RFA(K)=CRF*(A1-RRTP)
   15 IF(A2.GE.RRTP) RFA1(K)=CRF*(A2-RRTP)
      CKH=2.
      CKHB=2.
c
      A2K=0.2E6
c
      AAK=0.00175/DT
c
      if (ijkadv .eq. 1) then
         A2K=0.15E6 ! tao 5/1/01
c        A2K=0.2E6
c         aak=0.00175/dt	!tao-new
         aak=0.001250/dt  ! tao 5/1/01
      endif

      IF(ICPS_BM .EQ. 1) a2k=.1e6
      IF(ICPS_BM .EQ. 1) AAK=0.00125/DT
c
      kref=2
      zkref=250.    !m
      do k=2,kmax
        if(((z2(k)-z2(k-1))/100.).le.zkref) kref=k
      enddo
      print*,'kref is',kref
      A2I=AAK*DX*DX
      DO 20 K=1,KMAX
         Y3(K)=DZ*DZ/(AM(K)*AM(K))
         COEF(K)=CK*CK*Y3(K)*.5
         C1(K)=3.*980.*COEF(K)
         C2(K)=.5*CE/(CK*Y3(K))
         C3(K)=2./(3.*CK*CK*Y3(K))
ctao  (11-13)
        A4H(K)=AAK*(AM(kref)/AM(K))
c        A4H(K)=AAK
          IF(ICPS_BM .EQ. 1) A4H(K)=AAK*(AM(kref)/AM(K))
          A4K(K)=A4H(K)/CKHB
c
        TBSK(K)=2.0e6*(AM(kref)/AM(K)) ! tao 5/1/01
c       TBSK(K)=2.5e6*(AM(kref)/AM(K))
         IF(ICPS_BM .EQ. 1) TBSK(K)=2.5e6*(AM(kref)/AM(K))/10.
c
        BSKT(K)=A2K*(AM(kref)/AM(K))
          BSKT2(K)=2.*BSKT(K)
         BSKM(K)=BSKT(K)/CKH
           BSKM4(K)=4.*BSKM(K)
c
         BSIT(K)=A2I*(AM(kref)/AM(K))
           BSIT2(K)=2.*BSIT(K)
          BSIM(K)=BSIT(K)/CKH
            BSIM4(K)=4.*BSIM(K)
   20 CONTINUE
      if(isounding.eq.1)then 
        READ(5,202) (FD(K),TB(K),K=2,KMAX)
        READ(5,202) (FE(K),QB(K),K=2,KMAX)
        READ(5,250) (UB(K),K=2,KMAX)
        READ(5,250) (VB(K),K=2,KMAX)
        READ(5,202) (WB(K),K=2,KMAX)
        IF (IWBAR .EQ. 1) THEN
          DO ICY=1,itt
            READ(5,202) (WB_6H(K,ICY),K=2,KMAX)
          ENDDO
          IF (ISUI .EQ. 0) THEN
            DO ICY=1,itt
              READ(5,202) (Q1Z_6H(K,ICY),K=2,KMAX)
              READ(5,202) (Q2Z_6H(K,ICY),K=2,KMAX)
            ENDDO
          ENDIF
          DO K=2,KMAX
            QR(K)=WB(K)
            DO ICY=1,itt
              QR(K)=QR(K)+WB_6H(K,ICY)
            ENDDO
            QR(K)=QR(K)/(itt+1)
            WB(K)=QR(K)
          ENDDO
        ELSE
          DO ICY=1,itt
            READ(5,250) (UB_6H(K,ICY),K=2,KMAX)
            READ(5,250) (VB_6H(K,ICY),K=2,KMAX)
            READ(5,202) (WB_6H(K,ICY),K=2,KMAX)
          ENDDO
          do k=23,kmax
            ub(k)=ub(22)
            do icy=1,itt
              ub_6h(k,icy)=ub_6h(22,icy)
            enddo
          enddo
          DO ICY=1,itt
            READ(5,202) (Q1_6H(K,ICY),K=2,KMAX)
            READ(5,202) (Q1Z_6H(K,ICY),K=2,KMAX)
            READ(5,202) (Q2_6H(K,ICY),K=2,KMAX)
            READ(5,202) (Q2Z_6H(K,ICY),K=2,KMAX)
          ENDDO
          DO ICY=1,itt
            Q1Z_6H(1,ICY)=Q1Z_6H(2,ICY)
            Q2Z_6H(1,ICY)=Q2Z_6H(2,ICY)
          ENDDO
        ENDIF
      else
c
        call ainter(sst,ssp,tb,fd,qb,fe,ub,vb,wb,z1,z2,toptq,topuv)
c
        PSFC=SSP(1)*1000.
        do k=1,kmax
          if((z2(k)/100.).le.toptq) ktoptq=k
          if((z2(k)/100.).le.topuv) ktopuv=k
c          do i=1,itt
c            q1z_6h(k,i)=0.
c            q2z_6h(k,i)=0.
c          enddo
        enddo
       print*,'ktoptq, m, tophtq ',ktoptq,.01*z2(ktoptq),toptq
       print*,'ktopuv, m, tophuv ',ktopuv,.01*z2(ktopuv),topuv
      endif
c      rate=5.
c      rate=max(rate,tb(kmax)-tb(ktoptq))
c      rate=min(rate,20.)/(real(kmax-ktoptq)+1.e-10)
c      do k=ktoptq+1,kmax
c        tb(k)=rate+tb(k-1)
c        fd(k)=rate+fd(k-1)
c        qb(k)=0.
c        fe(k)=0.
c        wb(k)=0.
c      enddo
c
cc
ccc

cc
c
      IF (ILIF .EQ. 1) THEN
c       dili is used to read v,tvad,qvad
        r91=1.0
        r92=0.1
        z91=z1(18)
        z92=z1(26)
        a99=(r91-r92)/(z91-z92)
        b99=(r91*z92-r92*z91)/(z92-z91)
        do k=1,9
          z99=z1(17+k)
          r99=a99*z99+b99
         if(k.eq.1) ry1=r99
         if(k.eq.2) ry2=r99
         if(k.eq.3) ry3=r99
         if(k.eq.4) ry4=r99
         if(k.eq.5) ry5=r99
         if(k.eq.6) ry6=r99
         if(k.eq.7) ry7=r99
         if(k.eq.8) ry8=r99
         if(k.eq.9) ry9=r99
       enddo
       write(6,*) 'test',ry1,ry2,ry3,ry4,ry5,ry6,ry7,ry8,ry9
        wb(18)=ry1*wb(18)
	wb(19)=ry2*wb(19)
	wb(20)=ry3*wb(20)
	wb(21)=ry4*wb(21)
	wb(22)=ry5*wb(22)
	wb(23)=ry6*wb(23)
	wb(24)=ry7*wb(24)
	wb(25)=ry8*wb(25)
	wb(26)=ry9*wb(26)
	wb(27)=ry9*wb(27)
	wb(28)=ry9*wb(28)
	wb(29)=ry9*wb(29)
	wb(30)=ry9*wb(30)
	wb(31)=ry9*wb(31)
	wb(32)=ry9*wb(32)
	wb(33)=ry9*wb(33)
        do icy=1,ITT
           wb_6h(18,icy)=ry1*wb_6h(18,icy)
           wb_6h(19,icy)=ry2*wb_6h(19,icy)
	   wb_6h(20,icy)=ry3*wb_6h(20,icy)
	   wb_6h(21,icy)=ry4*wb_6h(21,icy)
 	   wb_6h(22,icy)=ry5*wb_6h(22,icy)
	   wb_6h(23,icy)=ry6*wb_6h(23,icy)
	   wb_6h(24,icy)=ry7*wb_6h(24,icy)
	   wb_6h(25,icy)=ry8*wb_6h(25,icy)
	   wb_6h(26,icy)=ry9*wb_6h(26,icy)
	   wb_6h(27,icy)=ry9*wb_6h(27,icy)
	   wb_6h(28,icy)=ry9*wb_6h(28,icy)
	   wb_6h(29,icy)=ry9*wb_6h(29,icy)
	   wb_6h(30,icy)=ry9*wb_6h(30,icy)
	   wb_6h(31,icy)=ry9*wb_6h(31,icy)
	   wb_6h(32,icy)=ry9*wb_6h(32,icy)
	   wb_6h(33,icy)=ry9*wb_6h(33,icy)
        enddo
        do k=23,kmax
          ub(k)=ub(22)
        do icy=1,ITT
          ub_6h(k,icy)=ub_6h(22,icy)
        enddo
        enddo
        do icy=1,ITT
	   q1_6h(18,icy)=ry1*q1_6h(18,icy)
	   q1_6h(19,icy)=ry2*q1_6h(19,icy)
	   q1_6h(20,icy)=ry3*q1_6h(20,icy)
	   q1_6h(21,icy)=ry4*q1_6h(21,icy)
	   q1_6h(22,icy)=ry5*q1_6h(22,icy)
	   q1_6h(23,icy)=ry6*q1_6h(23,icy)
	   q1_6h(24,icy)=ry7*q1_6h(24,icy)
	   q1_6h(25,icy)=ry8*q1_6h(25,icy)
	   q1_6h(26,icy)=ry9*q1_6h(26,icy)
	   q1_6h(27,icy)=ry9*q1_6h(27,icy)
	   q1_6h(28,icy)=ry9*q1_6h(28,icy)
	   q1_6h(29,icy)=ry9*q1_6h(29,icy)
	   q1_6h(30,icy)=ry9*q1_6h(30,icy)
	   q1_6h(31,icy)=ry9*q1_6h(31,icy)
	   q1_6h(32,icy)=ry9*q1_6h(32,icy)
	   q1_6h(33,icy)=ry9*q1_6h(33,icy)
cxx    do k=2,kmax
cxx       q1_6h(k,icy)=q1_6h(k,icy)+273.16
cxx       q2_6h(k,icy)=q2_6h(k,icy)*0.001
cxx    enddo
       enddo
       wb(1)=0.
       wb(2)=0.
       wb(kmax)=0.
       do icy=1,ITT
          ub_6h(   1,icy)=ub_6h(   2,icy)
          vb_6h(   1,icy)=vb_6h(   2,icy)
          q1_6h(   1,icy)=q1_6h(   2,icy)
          q2_6h(   1,icy)=q2_6h(   2,icy)
          wb_6h(   1,icy)=0.
          wb_6h(   2,icy)=0.
          wb_6h(kmax,icy)=0.
       end do
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,210)
      DO 32 K1=2,KMAX
        K=KMAX+2-K1
   32 WRITE(6,2031) K,TB(K),QB(K),FD(K),FE(K),QR(K),UB(K),VB(K),WB(K)
      WRITE(6,204)
      DO 34 K1=2,KMAX
        K=KMAX+2-K1
   34 WRITE(6,2441) K,Z1(K)/100.,Z2(K)/100.,AM(K),AM1(K),COEF(K)
     *,              C1(K),C2(K),C3(K),TBSK(K),BSKM(K),RFA(K)
      DO 40 K=2,KMAX
        FD(K)=FD(K)+273.16
        FE(K)=FE(K)*1.E-3
        TB(K)=TB(K)+273.16
   40   QB(K)=QB(K)*1.E-3
      A1=980.*DZ/RA
      P0(2)=LOG(PSFC)
      DO 45 K=3,KMAX
   45   P0(K)=P0(K-1)-A1/(TB(K-1)*(1.+.61*QB(K-1))*AM(K-1))

      RHO1(KMAX)=EXP(P0(KMAX))/(RA*FD(KMAX)*(1.+.61*FE(KMAX)))
      DO 46 K=2,KLES
   46   DZ0(K)=Z2(K+1)-Z2(K)
      DZ0(KMAX)=DZ0(KLES)
      DO 47 K=2,KMAX
   47   P00(K)=EXP(P0(K))
      DO 50 K=2,KLES
        RHO1(K)=EXP(P0(K))/(RA*FD(K)*(1.+.61*FE(K)))
        P0(K)=EXP((P0(K)+P0(K+1))*.5)
        PI(K)=(P0(K)/1000.E3)**.286
        F0(K)=AL/(PI(K)*CP)
        Y1(K)=TB(K)
        TB(K)=TB(K)/PI(K)
        QVS(K)=3.799E3/P0(K)*EXP(17.26939*(Y1(K)-273.16)/(Y1(K)-35.86))
        Y3(K)=3.799E3/P0(K)*EXP(21.87456-5807.695/(Y1(K)-7.66))
c        if(isounding.eq.1)then
c          if (k.eq.26) qb(k)=0.3*y3(k)
c          if (k.eq.27) qb(k)=0.225*y3(k)
c          if (k.eq.28) qb(k)=0.15*y3(k)
c          if (k.eq.29) qb(k)=0.015*y3(k)
c          if (k.eq.30) qb(k)=0.005*y3(k)
c          if (k.eq.31) qb(k)=0.001*y3(k)
c          if (k.ge.32) qb(k)=0.0   
c        endif
        Y4(K)=QB(K)/Y3(K)
        HU1(K)=QB(K)/QVS(K)
c tao add 5/01/01, try to resolve over-cooling problem
        if (k .ge. 23 .and. y4(k) .ge. 0.50) qb(k)=0.5*y3(k)
          if (k.eq.24) qb(k)=0.45*y3(k)
          if (k.eq.25) qb(k)=0.40*y3(k)
          if (k.eq.26) qb(k)=0.35*y3(k)
          if (k.eq.27) qb(k)=0.30*y3(k)
          if (k.eq.28) qb(k)=0.25*y3(k)
          if (k.eq.29) qb(k)=0.15*y3(k)
          if (k.eq.30) qb(k)=0.10*y3(k)
          if (k.eq.31) qb(k)=0.075*y3(k)
          if (k.ge.32) qb(k)=0.025*y3(k)  
          if (k.ge.33) qb(k)=0.0   
        Y4(K)=QB(K)/Y3(K)
        HU1(K)=QB(K)/QVS(K)

   50 RHO(K)=P0(K)/(RA*Y1(K)*(1.+.61*QB(K)))
      DO K=3,KLES-1
        IF(TB(K).LT.TB(K-1)) TB(K)=.5*(TB(K-1)+TB(K+1))
        IF(TB(K).LT.TB(K-1)) TB(K)=TB(K-1)-TB(K)+TB(K-1)
      ENDDO
      IF(TB(KLES).LT.TB(KLES-1))TB(KLES)=TB(KLES-1)*2.-TB(KLES-2)
        rqm=0.
        do k=2,kles
        rqm=rqm+(.5*(rho1(k)+rho(k))*.5*(qb(k)+fe(k)))*(hgt1(k)-hgt2(k))
        rqm=rqm+(.5*(rho1(k+1)+rho(k))*.5*(qb(k)+fe(k+1)))
     *     *(hgt2(k+1)-hgt1(k))
        enddo 
        print*,'MODEL QV INTEGRATION=',rqm
        print*

        if (ilif .eq. 0) then

         zpresst=75.*1000.
         zpressu=75.*1000.
        do k=2,kles
           if(p0(k) .ge. zpresst) toptq=k
           if(p0(k) .ge. zpressu) topuv=k
        enddo
        print*,'top at or just below 75mb =',toptq
        print*,'top at or just below 75mb =',topuv
          do k=ktopuv+1,kmax
            ub(k)=ub(ktopuv)
            vb(k)=vb(ktopuv)
          enddo
         do k=ktoptq+1,kmax
           wb(k)=0.
         enddo

cccshie  3/30/01 by dan: ADJUST THE VALUES NEAR THE TOP

         do k=2,kles
           p0mb=p0(k)*.001 
c          if (p0mb.lt.350.) then
           if (p0mb.lt.450.) then  ! cccshie by tao 5/3/01
              k350=k
              goto 4321
           endif
         enddo
4321     continue

         do k=2,kles
           p0mb=p0(k)*.001 
c          if (p0mb.lt.100.) then
           if (p0mb.lt.150.) then ! cccshie by tao 5/4/01
              k100=k
              goto 4323
           endif
         enddo
4323     continue

cccshie by tao 5/4/01
         do k=2,kles
           p0mb=p0(k)*.001
           if (p0mb.lt.250.) then
              k250=k
              goto 4325
           endif
         enddo
4325     continue

         do icy=1,itt
            do k=ktopuv+1,kmax
              ub_6h(k,icy)=ub_6h(ktopuv,icy)
              vb_6h(k,icy)=vb_6h(ktopuv,icy)
            enddo
            do k=ktoptq+1,kmax
              wb(k)=0.
              wb_6h(k,icy)=0.
c              q1_6h(k,icy)=0.   ! dan 3/30/01
c              q2_6h(k,icy)=0.   ! dan 3/30/01
c              q1z_6h(k,icy)=0.  ! dan 3/30/01
c              q2z_6h(k,icy)=0.  ! dan 3/30/01
            enddo
            wb_6h(2,icy)=0.
            wb_6h(1,icy)=0.
         
cccshie  3/30/01 by dan: ADJUST THE VALUES NEAR THE TOP
         do k=k350,k100
c           divfac=(1.-real(k-k350)/real(k100-k350))
            divfac=(1.-real(k-k350)/real(k100-k350))**2. ! cccshie by tao 5/3/01
            q1_6h(k,icy)=q1_6h(k350,icy)*divfac
c           q2_6h(k,icy)=q2_6h(k350,icy)*divfac    ! cccshie by tao 5/4/01
            q1z_6h(k,icy)=q1z_6h(k350,icy)*divfac
c           q2z_6h(k,icy)=q2z_6h(k350,icy)*divfac  ! cccshie by tao 5/4/01
         enddo
cccshie 5/4/01 shie modified following tao's idea: to drier top layers sounding
         do k=k350,k250
           divfacv=(1.-real(k-k350)/real(k250-k350))**2.
            q2_6h(k,icy)=q2_6h(k350,icy)*divfacv
            q2z_6h(k,icy)=q2z_6h(k350,icy)*divfacv
         enddo
         do k=k250,k100
            q2_6h(k,icy)=0.
            q2z_6h(k,icy)=0.
         enddo

         do k=k100+1,kmax
            q1_6h(k,icy)=0.
            q2_6h(k,icy)=0.
            q1z_6h(k,icy)=0.
            q2z_6h(k,icy)=0.
         enddo

       enddo
      endif

cc
cc
      IF (LIPPS.EQ.0)THEN
        DO K=3,KLES
          A1=RHO1(K)*.5
          FD(K)=(TB(K)-TB(K-1))*A1
          FE(K)=(QB(K)-QB(K-1))*A1
          A2=AM1(K)*(TB(K)+TB(K-1))*A1
          BA(K)=AM(K)*A2/(RHO(K)*TB(K))
          BB(K-1)=AM(K-1)*A2/(RHO(K-1)*TB(K-1))
        ENDDO
      ELSE
        DO K=3,KLES
          A1=RHO1(K)*.5
          FD(K)=(TB(K)-TB(K-1))*A1
          FE(K)=(QB(K)-QB(K-1))*A1
          A2=AM1(K)*RHO1(K)
          BA(K)=AM(K)*A2/RHO(K)
          BB(K-1)=AM(K-1)*A2/RHO(K-1)
        ENDDO
      ENDIF
      BA(2)=AM(2)*AM1(2)*RHO1(2)/RHO(2)
      BB(KLES)=AM(KLES)*AM1(KMAX)*RHO1(KMAX)/RHO(KLES)
c      if(iuvbar.eq.0)then 
c        A2=0.
c        A3=0.
c        A4=0.
c        DO 6666 K=2,KLES
c          A1=RHO(K)*DZ/AM(K)
c          A2=A2+A1*UB(K)
c          A3=A3+A1*VB(K)
c          A4=A4+A1
c 6666   CONTINUE
c        DO 6669 K=2,KLES
c 6669     VB(K)=A3/A4
c      endif
      do k=itt,1,-1
        if(isounding.eq.1)then
          TSL=(TSFC+273.16)/pi(2)
          QSL=3.799052E3/p0(2)*EXP(17.26939*TSFC/(TSFC+237.3))*.975
        else 
          TSFC=SST(k)
          PSFC=SSP(k)*1000.
          pisfc=(psfc/1000.e3)**.286
          TSL=(TSFC+273.16)/pisfc
          QSL=3.799052E3/psfc*EXP(17.26939*TSFC/(TSFC+237.3))
        endif
          print*,'psfc',psfc,psfc*.001
        write(6,12121) k,tsl,qsl*1000.,psfc*.001,tsfc+273.16
        do i=1,imax
          THS(i,k)=TSL
          QS(i,k)=QSL
          TS(i,k)=tsfc+273.16
          PSS(i,k)=psfc
        enddo
      enddo
12121 format(2x,3hit=,i4,4x,5htsfc=,f12.5,4x,5hqsfc=,f12.5,4x,
     1          5hpsfc=,f12.5,4x,7htsfc c=,f12.5)

      TB(1)=TB(2)
      QB(1)=QB(2)
      UB(1)=UB(2)
      VB(1)=VB(2)
      WB(1)=0.
      wb(kmax)=0.
      PI(1)=PI(2)
      HU1(1)=HU1(2)
      FD(2)=0.
      FE(2)=0.
      FD(1)=FD(2)
      FE(1)=FE(2)
      F0(1)=F0(2)
      RHO1(1)=RHO1(2)
      RHO(1)=RHO(2)
      P0(1)=P0(2)
      TB(KMAX)=TB(KLES)*2.-TB(KL2)
      QB(KMAX)=QB(KLES)*2.-QB(KL2)
      UB(KMAX)=UB(KLES)
      VB(KMAX)=VB(KLES)
      PI(KMAX)=PI(KLES)
      HU1(KMAX)=0.
      F0(KMAX)=2.*F0(KLES)-F0(KL2)
      FD(KMAX)=FD(KLES)
      FE(KMAX)=FE(KLES)
      RHO(KMAX)=RHO(KLES)
      P0(KMAX)=P0(KLES)
      DO 70 K=1,KMAX
        wb(k)=scale*wb(k)/(rho1(k)*980.)
        if(k.le.2) wb(k)=0.
        if(k.eq.kmax) wb(k)=0.
           if (p0(k) .le. 75.e3) then
             if (wb(k) .lt. 0.0) wb(k)=0.
           endif
         do icy=1,itt
          if (k.le.2) wb_6h(k,icy)=0.
          if (k.eq.kmax) wb_6h(k,icy)=0.
           wb_6h(k,icy)=scale*wb_6h(k,icy)/(rho1(k)*980.)
             if (p0(k) .le. 75.e3) then
               if (wb_6h(k,icy) .lt. 0.0) wb_6h(k,icy)=0.
             endif
         enddo
        RRHO(K)=1./RHO(K)
        RRHO1(K)=1./RHO1(K)
   70 CONTINUE
      IF(IRS.EQ.1) GO TO 1111
      DO 150 K=1,KMAX
        TA(K)=TB(K)
        QA(K)=QB(K)
        TA1(K)=TB(K)
        QA1(K)=QB(K)
        UB1(K)=UB(K)
        VB1(K)=VB(K)
      DO 150 I=1,IMAX
        UU1(I,K)=UB(K)
        VV1(I,K)=VB(K)
        ww1(i,k)=0.
        U(I,K)=UB(K)
        V(I,K)=VB(K)
        w(i,k)=0.
        umd(I,K)=UB(K)
        vmd(i,k)=vb(k)
        wmd(I,K)=0.
  150 CONTINUE
 1111 WRITE(6,200)
      DO 33 K1=2,KLES
      K=KLES+2-K1
      KM=K-1
   33 WRITE(6,2032)KM,P0(K)/1000.,PI(K),TB(K),QB(K)*1000.,QVS(K)*1000.
     *,     HU1(K),BA(K)*1000.,BB(K)*1000.,Y3(K)*1000.,Y4(K)
      WRITE(6,201)
      DO 35 K1=2,KLES
      K=KLES+2-K1
      KM=K-1
   35 WRITE(6,2033)KM,RHO(K)*1000.,RHO1(K)*1000.,F0(K),FD(K),FE(K),
     1             UB(K)/100.,VB(K)/100.,WB(K),BSIM(K)
      do icy=1,itt
        print*
        print*,'itt =',icy
        print*
        write(6,241)
        do 36 k1=2,kles
          k=kles+2-k1
          km=k-1
   36     write(6,203) km,P0(K)/1000.,ub_6h(k,icy),vb_6h(k,icy),
     1                   wb_6h(k,icy),Q1_6H(k,icy),Q2_6H(k,icy)*1000.,
     2                   t_adjust(k,icy),q_adjust(k,icy)*1000.
      enddo
c      stop
      RETURN
  200 FORMAT(//,' LEVEL   PRES    PI      TB       QB     QSW     HUM      
     *   BA      BB     QSI      HUI')
  201 FORMAT(//,'LEVEL  RHO   RHO1    F0        FD       FE       UB   
     1     VB    WB      BSIM')
  202 FORMAT(16F5.2)
  203 FORMAT(I4,10F12.5,/)
 2035 FORMAT(I4,2F10.5,4F12.3)
 2031 FORMAT(1X,I4,2(F8.2,F7.2),F8.2,2F10.2,F8.2)
 2032 FORMAT(I4,F9.2,F7.3,F9.2,7F8.2)
 2033 FORMAT(I4,2F7.3,F9.2,E11.3,E12.3,2F7.2,F7.3,E11.3)
  204 FORMAT(//,'LEVEL  Z      Z1     AM    AM1    COEF       C1       
     *C2       C3        TBSK       BSKM')
 2441 FORMAT(I3,2F8.1,2F6.3,5E10.3,2E9.2,E10.3)
  210 FORMAT(//,'  LEVEL  TB      QB     FD     FE1      QR       UB   
     1     VB        WB')
  241 format(//,'level          p           u           v           w 
     1       q1          q2        tdt       qdt')
  250 format(16f5.0)
      END
