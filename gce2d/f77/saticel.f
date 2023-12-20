
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SATICEL (n,id,fv)
C     (LIN)  COMPUTE ICE PHASE MICROPHYSICS AND SATURATION PROCESSES
      PARAMETER (NX=514,NZ=43,NT=2880,ITT=244)
      PARAMETER (NT5=5*NT,NB=NX*NZ-29*NX,NB1=NX*NZ-13*NX)

      integer n,id
      real    fv(nz)

      common/timestat/ ndt_stat,itime_ave,mmave
      integer itoga,ISFC,ICE,ICE2
      common/itoga/ itoga,ISFC,ICE,ICE2
      common/iceopt/ ice913,ilif
      common/ilin/ lin
      common/iice/ new_ice_sat
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bsize/ it1,it2,kt1,kt2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141

      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC
      COMMON/BSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,BND1,RN2,BND2,
     1   RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,RN101,RN10A,RN11,RN11A,
     2   RN12,RN12A(31),RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN17,
     3   RN17A,RN17B,RN17C,RN18,RN18A,RN19,RN19A,RN19B,RN20,RN20A,RN20B,
     4   BND3,RN21,RN22,RN23,RN23A,RN23B,RN25,RN25A(31),RN30A,RN30B,
     5   RN30C,RN31,BETA,RN32

      COMMON/B1T/ DPT(NX,NZ)
      COMMON/B1Q/ DQV(NX,NZ)
      COMMON/B1C/ QCL(NX,NZ)
      COMMON/B1R/ QRN(NX,NZ)
      COMMON/B2T/ DPT1(NX,NZ)
      COMMON/B2Q/ DQV1(NX,NZ)
      COMMON/B2C/ QCL1(NX,NZ)
      COMMON/B2R/ QRN1(NX,NZ)
      COMMON/B1I/ QCI(NX,NZ)
      COMMON/B1S/ QCS(NX,NZ)
      COMMON/B1G/ QCG(NX,NZ)
      COMMON/B2I/ QCI1(NX,NZ)
      COMMON/B2S/ QCS1(NX,NZ)
      COMMON/B2G/ QCG1(NX,NZ)
      COMMON/B2W/ WW1(NX,NZ)
      COMMON/SLWAVE/ RSW(NX,NZ),RLW(NX,NZ)

      COMMON/BSAT/ WGACR(NX),DEP(NX),RGMP(NX),
     1   DD(NX),DD1(NX),QVS(NX),DM(NX),RSUB1(NX),COL(NX),CND(NX),RQ(NX),
     2   ERN(NX),SCV(NX),TCA(NX),DWV(NX),ZR(NX),VR(NX),ZS(NX),VS(NX),
     3   ZG(NX),VG(NX),EGS(NX),ESI(NX),QSI(NX),SSI(NX),QSW(NX),SSW(NX),
     4   PIHOM(NX),PIDW(NX),DDC(NB)

      COMMON/BSAT1/ QSACW(NX),PRACI(NX),PIACR(NX),PRAUT(NX),
     1   PRACW(NX),PSFW(NX),PSFI(NX),DGACS(NX),DGACW(NX),DGACI(NX),
     2   DGACR(NX),PGACS(NX),WGACS(NX),QGACW(NX),WGACI(NX),QGACR(NX),
     3   PGWET(NX),PGAUT(NX),PRACS(NX),PSACR(NX),QSACR(NX),PGFR(NX),
     4   PSMLT(NX),PGMLT(NX),PSDEP(NX),PSSUB(NX),PGSUB(NX),PINT(NX),
     5   PIDEP(NX),DDA(NB)


      COMMON/BADV/ PT(NX),QV(NX),QC(NX),QR(NX),QI(NX),QS(NX),QG(NX),
     1   XX0(NX),XX00(NX),PIMLT(NX),PSAUT(NX),PSACI(NX),PSACW(NX),
     2   DDB(NB1)

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),tair(nx),tairc(nx),
     $   pr(nx),ps(nx),pg(nx),prn(nx),psn(nx),dlt1(nx),dlt2(nx),
     $   dlt3(nx),rtair(nx)
      common/ba/ y1,y2,y3,y4,y5,tair,tairc,pr,ps,pg,prn,psn,dlt1,dlt2,
     $           dlt3,rtair
c
c
      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     $   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    srro(nz),qrro(nz),sqc(nz),sqr(nz),sqi(nz),sqs(nz),
     $   sqg(nz),stqc(nz),stqr(nz),stqi(nz),stqs(nz),stqg(nz),ttt(nt5)
      common/brh1/ srro,qrro,sqc,sqr,sqi,sqs,
     $             sqg,stqc,stqr,stqi,stqs,stqg,ttt

      COMMON/B7/ SQAAQ(NZ),SQTT(NZ),SQVV(NZ),SQAAT(NZ),SQAK(NZ)
      real    y0(nx),ths(nx,itt),qsTT(nx,itt),tsTT(nx,itt),pss(nx,itt)
      common/bls/ y0,ths,qsTT,tsTT,pss

      integer it(nx),ics(nx,4),ibz(nx,4),iv(nt)
      common/bch/ it,iv,ics,ibz
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON/BSTS/ THOM(NZ,4,7),TDW(NZ,4,7),TMLT(NZ,4,7),
     1 SAUT(NZ,4,7),SACI(NZ,4,7),SACW(NZ,4,7),RACI(NZ,4,7),
     2 TACR(NZ,4,7),RAUT(NZ,4,7),RACW(NZ,4,7),SFW(NZ,4,7),
     3 SFI(NZ,4,7),GACS(NZ,4,7),GACW(NZ,4,7),
     4 GACI(NZ,4,7),GACR(NZ,4,7),GWET(NZ,4,7),GAUT(NZ,4,7),
     5 RACS(NZ,4,7),SACR(NZ,4,7),GFR(NZ,4,7),SMLT(NZ,4,7),
     6 GMLT(NZ,4,7),SDEP(NZ,4,7),SSUB(NZ,4,7),GSUB(NZ,4,7),
     7 PERN(NZ,4,7),D3RI(NZ,4,7),D3IR(NZ,4,7),D2SR(NZ,4,7),
     8 D2RS(NZ,4,7),GDRY(NZ,4,7),COC(NZ,4,7),COE(NZ,4,7),
     9 SMF0(NZ,4,7),QC0(NZ,4,7),QR0(NZ,4,7),QI0(NZ,4,7),
     1 QS0(NZ,4,7),QG0(NZ,4,7),SQC0(NZ,4,7),SQR0(NZ,4,7),
     2 SQI0(NZ,4,7),SQS0(NZ,4,7),SQG0(NZ,4,7),ERNS(NZ,4,7),
     3 WGRS(NZ,4,7),QSWS(NZ,4,7),TBT(NZ,4),QBT(NZ,4)


      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7)

      COMMON/BSTSLIN/ qv0l(nz,4,7),tt0l(nz,4,7),sgpt(nz,4,7),
     1  tsqq(nz,4,7),stt0l(nz,4,7),sqv0l(nz,4,7),ssgpt(nz,4,7),
     2  tsqq1(nz,4,7),q1tl(nz,4,7),sqqdtl(nz,4,7),sqdtl(nz,4,7),
     3  stdtl(nz,4,7)

      real   S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     $       S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     $       S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     $       SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     $       SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     $       SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCS/ S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,
     $            S19,S20,S21,SCNT,SN9,SN10,SN11,SN12,SN13,SN14,
     $            SN15,SN16,SN17,SN18,SN19,SN20,SN21,SNCNT,SCU1,SED1

      real    SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     $ SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     $ SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     $ SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     $ SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     $ SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)
      COMMON/BCSS/ SS9,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,
     $         SS18,SS19,SS20,SS21,SSCNT,SSN9,SSN10,SSN11,SSN12,SSN13,
     $         SSN14,SSN15,SSN16,SSN17,SSN18,SSN19,SSN20,SSN21,SSNCNT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    rby(7)
      common/bch1/ rby

      DIMENSION AA1(31),AA2(31),dbz(nx)
      DATA AA1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,
     1   .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,
     2   .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,
     3   .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,
     4   .6513E-6,.5956E-6,.5333E-6,.4834E-6/
      DATA AA2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
      save
C     ******   THREE CLASSES OF ICE-PHASE   (LIN ET AL, 1983)  *********
       RT0=1./(T0-T00)
       D22T=D2T
       IF(IJKADV .EQ. 0) THEN
         D2T=D2T
       ELSE
         D2T=DT
       ENDIF
c
        UCOR=3071.29/TNW**.75                                            
        UCOS=687.97*ROQS**.25/TNS**.75
        UCOG=687.97*ROQG**.25/TNG**.75
        UWET=4.464**.95
       cmin=1.e-20
       cmin1=1.e-12
c
       A0=ndt_stat*.5*RIL2
      DO 10 I=1,imax
   10  IT(I)=1
       FT=DT/D2T
       RFT=RIL2*FT
        BW3=BW+3.
        BS3=BS+3.
        BG3=BG+3.
        BWH5=2.5+BWH
        BSH5=2.5+BSH
        BGH5=2.5+BGH
        BW6=BW+6.
        BETAH=.5*BETA
       R10T=RN10*D2T
       R11AT=RN11A*D2T
       R19BT=RN19B*D2T
       R20T=-RN20*D2T
       R23T=-RN23*D2T
       R25A=RN25
CC    ******************************************************************
      DO 1000 K=2,KLES
c        if (ijkadv .eq. 1) then
c          tb0=ta(k)
c          qb0=qa(k)
c        else
          tb0=ta1(k)
          qb0=qa1(k)
c        endif
c       P00=P0(K)
       RP0=3.799052E3/P0(K)
       PI0=PI(K)
       PIR=1./(PI(K))
       PR0=1./P0(K)
       R00=RHO(K)
       RR0=RRHO(K)
       RRS=SRRO(K)
       RRQ=QRRO(K)
       FV0=FV(K)
       FVS=SQRT(FV(K))
         ZRR=1.E5*ZRC*RRQ
         ZSR=1.E5*ZSC*RRQ
         ZGR=1.E5*ZGC*RRQ
         CP409=C409*PI0
         CV409=C409*AVC
         CP580=C580*PI0
         CS580=C580*ASC
         ALVR=R00*ALV
         AFCP=AFC*PIR
         AVCP=AVC*PIR
         ASCP=ASC*PIR
         VRCF=VRC*FV0
         VSCF=VSC*FV0

         VGCR=VGC*RRS
         VGCF=VGC*FV0

         DWVP=C879*PR0
         R3F=RN3*FV0
         R4F=RN4*FV0
         R5F=RN5*FV0
         R6F=RN6*FV0
         R7R=RN7*RR0
         R8R=RN8*RR0
         R9R=RN9*RR0
         R101F=RN101*FVS
         R10AR=RN10A*R00
         R11RT=RN11*RR0*D2T
         R12R=RN12*R00
         R14R=RN14*RRS
         R15R=RN15*RRS
         R15AR=RN15A*RRS
         R16R=RN16*RR0
         R17R=RN17*RR0
         R17AQ=RN17A*RRQ
         R18R=RN18*RR0
         R19RT=RN19*RR0*D2T
         R19AQ=RN19A*RRQ
         R20BQ=RN20B*RRQ
         R22F=RN22*FV0
         R23AF=RN23A*FVS
         R23BR=RN23B*R00
         R25RT=RN25*RR0*D2T
         R31R=RN31*RR0
         R32RT=RN32*D2T*RRS
       DO 100 I=2,ILES
        PT(I)=DPT(I,K)
        QV(I)=DQV(I,K)
        QC(I)=QCL(I,K)
        QR(I)=QRN(I,K)
        QI(I)=QCI(I,K)
        QS(I)=QCS(I,K)
        QG(I)=QCG(I,K)
c        IF (QV(I)+QB0 .LE. 0.) QV(I)=-QB0
        IF (QC(I) .LE. cmin) QC(I)=0.
        IF (QR(I) .LE. cmin) QR(I)=0.
        IF (QI(I) .LE. cmin) QI(I)=0.
        IF (QS(I) .LE. cmin) QS(I)=0.
        IF (QG(I) .LE. cmin) QG(I)=0.
        TAIR(I)=(PT(I)+TB0)*PI0
        TAIRC(I)=TAIR(I)-T0
C     ***   COMPUTE ZR,ZS,ZG,VR,VS,VG      *****************************
          ZR(I)=ZRR
          ZS(I)=ZSR
          ZG(I)=ZGR
          VR(I)=0.0
          VS(I)=0.0
          VG(I)=0.0
           IF (QR(I) .GT. cmin) THEN
             DD(I)=R00*QR(I)
              Y1(I)=DD(I)**.25
            ZR(I)=ZRC/Y1(I)
            VR(I)=min(max(VRCF*DD(I)**BWQ, 0.0), 1500.)
           ENDIF
           IF (QS(I) .GT. cmin) THEN
             DD(I)=R00*QS(I)
              Y1(I)=DD(I)**.25
            ZS(I)=ZSC/Y1(I)
            VS(I)=min(max(VSCF*DD(I)**BSQ, 0.0), 400.)
           ENDIF
           IF (QG(I) .GT. cmin) THEN
             DD(I)=R00*QG(I)
              Y1(I)=DD(I)**.25
            ZG(I)=ZGC/Y1(I)
            IF (LIN .EQ. 1) THEN
             vg(i)=min(max(vgcr*dd(i)**bgq, 0.), 1500.)
            ELSE
             vg(i)=min(max(vgcf*dd(i)**bgq, 0.), 800.)
            ENDIF
           ENDIF
         IF (QR(I) .LE. cmin1) VR(I)=0.0
         IF (QS(I) .LE. cmin1) VS(I)=0.0
         IF (QG(I) .LE. cmin1) VG(I)=0.0
C     ******************************************************************
C     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
C     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
C     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
C     ***   Y2 : KINETIC VISCOSITY (V)
          Y1(I)=C149*TAIR(I)**1.5/(TAIR(I)+120.)
        DWV(I)=DWVP*TAIR(I)**1.81
        TCA(I)=C141*Y1(I)
        SCV(I)=1./((RR0*Y1(I))**.1666667*DWV(I)**.3333333)
  100  CONTINUE
C*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
C*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
C*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
C*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
C*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
        DO 125 I=2,ILES
         PSAUT(I)=0.0
         PSACI(I)=0.0
         PRACI(I)=0.0
         PIACR(I)=0.0
         PSACW(I)=0.0
         QSACW(I)=0.0
          DD(I)=1./ZS(I)**BS3
         IF (TAIR(I).LT.T0) THEN
           ESI(I)=EXP(.025*TAIRC(I))
          PSAUT(I)=MAX(RN1*ESI(I)*(QI(I)-BND1) ,0.0)
          PSACI(I)=R3F*ESI(I)*QI(I)*DD(I)
          PSACW(I)=R4F*QC(I)*DD(I)
          PRACI(I)=R5F*QI(I)/ZR(I)**BW3
          PIACR(I)=R6F*QI(I)/ZR(I)**BW6
         ELSE
          QSACW(I)=R4F*QC(I)*DD(I)
         ENDIF
C* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
C* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
CM        PRAUT(I)=MAX((RD1*(QC(I)-BOUND)), 0.0)
CM        PRACW(I)=0.0
CM         IF (QR(I).GT.1.E-40) THEN
CM          PRACW(I)=RD2*QC(I)*QR(I)**.875
CM         ENDIF
          PRACW(I)=R22F*QC(I)/ZR(I)**BW3
          PRAUT(I)=0.0
            Y1(I)=QC(I)-BND3
           IF (Y1(I).GT.0.0) THEN
            PRAUT(I)=R00*Y1(I)*Y1(I)/(1.2E-4+RN21/Y1(I))
           ENDIF
C* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
C* 13 * PSFI : BERGERON PROCESSES FOR QS                          **13**
          PSFW(I)=0.0
          PSFI(I)=0.0
           IF(TAIR(I).LT.T0) THEN
             Y1(I)=MAX( MIN(TAIRC(I), -1.), -31.)
             IT(I)=INT(ABS(Y1(I)))
             Y1(I)=RN12A(IT(I))
             Y2(I)=RN12B(IT(I))
            PSFW(I)=MAX(D2T*Y1(I)*(Y2(I)+R12R*QC(I))*QI(I),0.0)
            PSFI(I)=RN13(IT(I))*QI(I)
           ENDIF
CTTT***** QG=QG+MIN(PGDRY,PGWET)
C*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
C* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
C* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
           EE1=1.
           EE2=0.09
          EGS(I)=EE1*EXP(EE2*TAIRC(I))
           IF (TAIR(I).GE.T0) EGS(I)=1.0
          Y1(I)=ABS(VG(I)-VS(I))
          Y2(I)=ZS(I)*ZG(I)
          Y3(I)=5./Y2(I)
          Y4(I)=.08*Y3(I)*Y3(I)
          Y5(I)=.05*Y3(I)*Y4(I)
          DD(I)=Y1(I)*(Y3(I)/ZS(I)**5+Y4(I)/ZS(I)**3+Y5(I)/ZS(I))
         PGACS(I)=R9R*EGS(I)*DD(I)
         DGACS(I)=PGACS(I)
         WGACS(I)=R9R*DD(I)
          Y1(I)=1./ZG(I)**BG3
         DGACW(I)=MAX(R14R*QC(I)*Y1(I), 0.0)
         QGACW(I)=DGACW(I)
          Y1(I)=ABS(VG(I)-VR(I))
          Y2(I)=ZR(I)*ZG(I)
          Y3(I)=5./Y2(I)
          Y4(I)=.08*Y3(I)*Y3(I)
          Y5(I)=.05*Y3(I)*Y4(I)
           DD(I)=R16R*Y1(I)*(Y3(I)/ZR(I)**5+Y4(I)/ZR(I)**3+Y5(I)/ZR(I))
         DGACR(I)=MAX(DD(I), 0.0)
         QGACR(I)=DGACR(I)
         IF (TAIR(I).GE.T0) THEN
          DGACS(I)=0.0
          WGACS(I)=0.0
          DGACW(I)=0.0
          DGACR(I)=0.0
         ELSE
          PGACS(I)=0.0
          QGACW(I)=0.0
          QGACR(I)=0.0
         ENDIF
C*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
C* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
C* 17 * PGWET : WET GROWTH OF QG                                  **17**
          DGACI(I)=0.0
          WGACI(I)=0.0
          PGWET(I)=0.0
          IF (TAIR(I).LT.T0 .AND. TAIR(IJ) .GE. 253.16) THEN
             Y1(I)=QI(I)/ZG(I)**BG3
            DGACI(I)=R15R*Y1(I)
            WGACI(I)=R15AR*Y1(I)
             Y1(I)=1./(ALF+RN17C*TAIRC(I))
             Y3(I)=.78/ZG(I)**2+R17AQ*SCV(I)/ZG(I)**BGH5
             Y4(I)=ALVR*DWV(I)*(RP0-(QV(I)+QB0))-TCA(I)*TAIRC(I)
            DD(I)=Y1(I)*(R17R*Y4(I)*Y3(I)
     1                      +(WGACI(I)+WGACS(I))*(ALF+RN17B*TAIRC(I)))
           PGWET(I)=MAX(DD(I), 0.0)
          ENDIF
  125   CONTINUE
C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
C********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
        DO 150 I=2,ILES

           Y1(I)=QC(I)/D2T
          PSACW(I)=MIN(Y1(I), PSACW(I))
          PRAUT(I)=MIN(Y1(I), PRAUT(I))
          PRACW(I)=MIN(Y1(I), PRACW(I))
          PSFW(I)= MIN(Y1(I), PSFW(I))
          DGACW(I)=MIN(Y1(I), DGACW(I))
          QSACW(I)=MIN(Y1(I), QSACW(I))
          QGACW(I)=MIN(Y1(I), QGACW(I))

        Y1(I)=(PSACW(I)+PRAUT(I)+PRACW(I)+PSFW(I)+DGACW(I)+QSACW(I)
     1         +QGACW(I))*D2T
       QC(I)=QC(I)-Y1(I)
        if (qc(i).lt.0.0) then
          a1=1.
           if(y1(i) .ne. 0.) a1=qc(i)/y1(i)+1.
         psacw(i)=psacw(i)*a1
         praut(i)=praut(i)*a1
         pracw(i)=pracw(i)*a1
         psfw(i)=psfw(i)*a1
         dgacw(i)=dgacw(i)*a1
         qsacw(i)=qsacw(i)*a1
         qgacw(i)=qgacw(i)*a1
         qc(i)=0.0
        endif



            Y1(I)=QI(I)/D2T
           PSAUT(I)=MIN(Y1(I), PSAUT(I))
           PSACI(I)=MIN(Y1(I), PSACI(I))
           PRACI(I)=MIN(Y1(I), PRACI(I))
           PSFI(I)= MIN(Y1(I), PSFI(I))
           DGACI(I)=MIN(Y1(I), DGACI(I))
           WGACI(I)=MIN(Y1(I), WGACI(I))

        Y1(I)=(PSAUT(I)+PSACI(I)+PRACI(I)+PSFI(I)+DGACI(I)+WGACI(I))*D2T
        QI(I)=QI(I)-Y1(I)

        if (qi(i).lt.0.0) then
          a2=1.
           if(y1(i) .ne. 0.) a2=qi(i)/y1(i)+1.
         psaut(i)=psaut(i)*a2
         psaci(i)=psaci(i)*a2
         praci(i)=praci(i)*a2
         psfi(i)=psfi(i)*a2
         dgaci(i)=dgaci(i)*a2
         wgaci(i)=wgaci(i)*a2
         qi(i)=0.0
        endif
c
C******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
         WGACR(I)=PGWET(I)-DGACW(I)-WGACI(I)-WGACS(I)
          Y2(I)=DGACW(I)+DGACI(I)+DGACR(I)+DGACS(I)
         IF (PGWET(I).GE.Y2(I)) THEN
          WGACR(I)=0.0
          WGACI(I)=0.0
          WGACS(I)=0.0
         ELSE
          DGACR(I)=0.0
          DGACI(I)=0.0
          DGACS(I)=0.0
         ENDIF
c******************************************************
c
         DLT3(I)=0.0
         DLT2(I)=0.0
          IF (TAIR(I).LT.T0) THEN
           IF (QR(I).LT.1.E-4) THEN
             DLT3(I)=1.0
             DLT2(I)=1.0
           ENDIF
           IF (QS(I).GE.1.E-4) THEN
             DLT2(I)=0.0
           ENDIF
          ENDIF
         PR(I)=(QSACW(I)+PRAUT(I)+PRACW(I)+QGACW(I))*D2T
         PS(I)=(PSAUT(I)+PSACI(I)+PSACW(I)+PSFW(I)+PSFI(I)
     1          +DLT3(I)*PRACI(I))*D2T
         PG(I)=((1.-DLT3(I))*PRACI(I)+DGACI(I)+WGACI(I)+DGACW(I))*D2T
  150   CONTINUE
C*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
C*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
        DO 175 I=2,ILES
          Y1(I)=ABS(VR(I)-VS(I))
          Y2(I)=ZR(I)*ZS(I)
          Y3(I)=5./Y2(I)
          Y4(I)=.08*Y3(I)*Y3(I)
          Y5(I)=.05*Y3(I)*Y4(I)
         PRACS(I)=R7R*Y1(I)*(Y3(I)/ZS(I)**5+Y4(I)/ZS(I)**3+Y5(I)/ZS(I))
         PSACR(I)=R8R*Y1(I)*(Y3(I)/ZR(I)**5+Y4(I)/ZR(I)**3+Y5(I)/ZR(I))
         QSACR(I)=PSACR(I)
         IF (TAIR(I).GE.T0) THEN
          PRACS(I)=0.0
          PSACR(I)=0.0
         ELSE
          QSACR(I)=0.0
         ENDIF
C*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
C* 18 * PGFR : FREEZING OF QR TO QG                               **18**
         PGAUT(I)=0.0
         PGFR(I)=0.0
          IF (TAIR(I) .LT. T0) THEN
            Y1(I)=EXP(.09*TAIRC(I))
            IF (LIN .EQ. 1) PGAUT(IJ)=MAX(RN2*Y1(IJ)*(QS(IJ)-BND2), 0.0)
            Y2(I)=EXP(RN18A*(T0-TAIR(I)))
           PGFR(I)=MAX(R18R*(Y2(I)-1.)/ZR(I)**7., 0.0)
          ENDIF
  175   CONTINUE
C********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
C********   HANDLING THE NEGATIVE SNOW (QS)          *******************
        DO 200 I=2,ILES

          Y1(I)=QR(I)/D2T
         PIACR(I)=MIN(Y1(I), PIACR(I))
         DGACR(I)=MIN(Y1(I), DGACR(I))
         WGACR(I)=MIN(Y1(I), WGACR(I))
         PSACR(I)=MIN(Y1(I), PSACR(I))
         PGFR(I)= MIN(Y1(I), PGFR(I))

          Y1(I)=(PIACR(I)+DGACR(I)+WGACR(I)+PSACR(I)+PGFR(I))*D2T
         QR(I)=QR(I)+PR(I)-Y1(I)
         if (qr(i).lt.0.0) then
           a1=1.
            if(y1(i) .ne. 0.) a1=qr(i)/y1(i)+1.
          PIACR(I)=PIACR(I)*A1
          DGACR(I)=DGACR(I)*A1
          WGACR(I)=WGACR(I)*A1
          PGFR(I)=PGFR(I)*A1
          PSACR(I)=PSACR(I)*A1
          QR(I)=0.0
         ENDIF
          PRN(I)=D2T*((1.-DLT3(I))*PIACR(I)+DGACR(I)+WGACR(I)
     1                +(1.-DLT2(I))*PSACR(I)+PGFR(I))
          PS(I)=PS(I)+D2T*(DLT3(I)*PIACR(I)+DLT2(I)*PSACR(I))
           PRACS(I)=(1.-DLT2(I))*PRACS(I)

           Y1(I)=QS(I)/D2T
          PGACS(I)=MIN(Y1(I), PGACS(I))
          DGACS(I)=MIN(Y1(I), DGACS(I))
          WGACS(I)=MIN(Y1(I), WGACS(I))
          PGAUT(I)=MIN(Y1(I), PGAUT(I))
          PRACS(I)=MIN(Y1(I), PRACS(I))

          PSN(I)=D2T*(PGACS(I)+DGACS(I)+WGACS(I)+PGAUT(I)+PRACS(I))
         QS(I)=QS(I)+PS(I)-PSN(I)
          if (qs(i).lt.0.) then
            A2=1.
             IF (PSN(I) .NE. 0.0) a2=qs(i)/psn(i)+1.
           pgacs(i)=pgacs(i)*a2
           dgacs(i)=dgacs(i)*a2
           wgacs(i)=wgacs(i)*a2
           pgaut(i)=pgaut(i)*a2
           pracs(i)=pracs(i)*a2
c           psn(i)=psn(i)*a2
           qs(i)=0.0
          endif
c
          PSN(I)=D2T*(PGACS(I)+DGACS(I)+WGACS(I)+PGAUT(I)+PRACS(I))
           Y2(I)=D2T*(PSACW(I)+PSFW(I)+DGACW(I)+PIACR(I)+DGACR(I)
     1                +WGACR(I)+PSACR(I)+PGFR(I))
         PT(I)=PT(I)+AFCP*Y2(I)
         QG(I)=QG(I)+PG(I)+PRN(I)+PSN(I)
  200   CONTINUE
C* 11 * PSMLT : MELTING OF QS                                     **11**
C* 19 * PGMLT : MELTING OF QG TO QR                               **19**
        DO 225 I=2,ILES
        PSMLT(I)=0.0
        PGMLT(I)=0.0
           TAIR(I)=(PT(I)+TB0)*PI0
        IF (TAIR(I).GE.T0) THEN
           TAIRC(I)=TAIR(I)-T0
           Y1(I)=TCA(I)*TAIRC(I)-ALVR*DWV(I)*(RP0-(QV(I)+QB0))
           Y2(I)=.78/ZS(I)**2+R101F*SCV(I)/ZS(I)**BSH5
           DD(I)=R11RT*Y1(I)*Y2(I)+R11AT*TAIRC(I)*(QSACW(I)+QSACR(I))
         PSMLT(I)=MAX(0.0, MIN(DD(I), QS(I)))
           Y3(I)=.78/ZG(I)**2+R19AQ*SCV(I)/ZG(I)**BGH5
           DD1(I)=R19RT*Y1(I)*Y3(I)+R19BT*TAIRC(I)*(QGACW(I)+QGACR(I))
         PGMLT(I)=MAX(0.0, MIN(DD1(I), QG(I)))
C         IF(ICS5(I,3,1).EQ.1 .OR. ICS5(I,4,1).EQ.1) THEN
C         IF(ICS5(I,2,1).EQ.1) THEN
C           PSMLT(I)=0.0
C           PGMLT(I)=0.0
C         ENDIF
         PT(I)=PT(I)-AFCP*(PSMLT(I)+PGMLT(I))
         QR(I)=QR(I)+PSMLT(I)+PGMLT(I)
         QS(I)=QS(I)-PSMLT(I)
         QG(I)=QG(I)-PGMLT(I)
        ENDIF

C* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
C* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0 < T <= T00)     **25**
C* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
        IF (QC(I).LE.cmin) QC(I)=0.0
        IF (QI(I).LE.cmin) QI(I)=0.0
          TAIR(I)=(PT(I)+TB0)*PI0
         PIHOM(I)=CVMGP(0.,QC(I),TAIR(I)-T00)
         PIMLT(I)=CVMGP(QI(I),0.,TAIR(I)-T0)
         PIDW(I)=0.0
         IF (TAIR(I).LT.T0 .AND. TAIR(I).GT.T00) THEN
          TAIRC(I)=TAIR(I)-T0
            Y1(I)=MAX( MIN(TAIRC(I), -1.), -31.)
            IT(I)=INT(ABS(Y1(I)))
CLIN        Y2(I)=RN25A(IT(I))
CLIN        Y3(I)=EXP(.5*ABS(TAIRC(I)))
CLIN       PIDW(I)=MIN(R25RT*Y2(I)*Y3(I), QC(I))
C
           Y2(I)=AA1(IT(I))
           Y3(I)=AA2(IT(I))
           Y4(I)=EXP(.5*ABS(TAIRC(I)))
           DD(I)=(R00*QI(I)/(R25A*Y4(I)))**Y3(I)
          PIDW(I)=MIN(R25RT*Y2(I)*Y4(I)*DD(I), QC(I))
         ENDIF
          Y1(I)=PIHOM(I)-PIMLT(I)+PIDW(I)
        PT(I)=PT(I)+AFCP*Y1(I)
        QC(I)=QC(I)-Y1(I)
        QI(I)=QI(I)+Y1(I)
C* 31 * PINT  : INITIATION OF QI                                  **31**
C* 32 * PIDEP : DEPOSITION OF QI                                  **32**
        PINT(I)=0.0
         TAIR(I)=(PT(I)+TB0)*PI0
         IF (TAIR(I) .LT. T0) THEN
            if (qi(i) .le. cmin) qi(i)=0.
           TAIRC(I)=TAIR(I)-T0
           DD(I)=R31R*EXP(BETA*TAIRC(I))
            RTAIR(I)=1./(TAIR(I)-C76)
             Y2(I)=EXP(C218-C580*RTAIR(I))
            QSI(I)=RP0*Y2(I)
            ESI(I)=C610*Y2(I)
            SSI(I)=(QV(I)+QB0)/QSI(I)-1.
           DM(I)=MAX( (QV(I)+QB0-QSI(I)), 0.)
            RSUB1(I)=CS580*QSI(I)*RTAIR(I)*RTAIR(I)
          PINT(I)=max(MIN(DD(I),DM(I)), 0.)
             Y1(I)=1./TAIR(I)
             Y2(I)=EXP(BETAH*TAIRC(I))
             Y3(I)=SQRT(QI(I))
           DD(I)=Y1(I)*(RN30A*Y1(I)-RN30B)+RN30C*TAIR(I)/ESI(I)
          PIDEP(I)=MAX(R32RT*SSI(I)*Y2(I)*Y3(I)/DD(I), 0.0)
          PINT(I)=PINT(I)+PIDEP(I)
           DEP(I)=DM(I)/(1.+RSUB1(I))
           PINT(I)=MIN(PINT(I), DEP(I))
cc            if (PINT(i) .le. cmin) PINT(i)=0.
          PT(I)=PT(I)+ASCP*PINT(I)
          QV(I)=QV(I)-PINT(I)
          QI(I)=QI(I)+PINT(I)
         ENDIF
  225   CONTINUE
C*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
        DO 250 I=2,ILES
         TAIR(I)=(PT(I)+TB0)*PI0
        CND(I)=RT0*(TAIR(I)-T00)
        DEP(I)=RT0*(T0-TAIR(I))
          Y1(I)=1./(TAIR(I)-C358)
          Y2(I)=1./(TAIR(I)-C76)
         QSW(I)=RP0*EXP(C172-C409*Y1(I))
         QSI(I)=RP0*EXP(C218-C580*Y2(I))
          DD(I)=CP409*Y1(I)*Y1(I)
          DD1(I)=CP580*Y2(I)*Y2(I)
         IF (QC(I).LE.cmin) QC(I)=cmin
         IF (QI(I).LE.cmin) QI(I)=cmin
         IF (TAIR(I).GE.T0) THEN
          DEP(I)=0.0
          CND(I)=1.0
          QI(I)=0.0
         ENDIF
         IF (TAIR(I).LE.T00) THEN
          CND(I)=0.0
          DEP(I)=1.0
          QC(I)=0.0
         ENDIF
          Y5(I)=AVCP*CND(I)+ASCP*DEP(I)
           Y1(I)=QC(I)*QSW(I)/(QC(I)+QI(I))
           Y2(I)=QI(I)*QSI(I)/(QC(I)+QI(I))
          Y4(I)=DD(I)*Y1(I)+DD1(I)*Y2(I)
         QVS(I)=Y1(I)+Y2(I)
         DM(I)=QV(I)+QB0-QVS(I)
         RSUB1(I)=DM(I)/(1.+Y4(I)*Y5(I))
        CND(I)=CND(I)*RSUB1(I)
        DEP(I)=DEP(I)*RSUB1(I)
         IF (QC(I).LE.cmin) QC(I)=0.
         IF (QI(I).LE.cmin) QI(I)=0.
CC    ******   CONDENSATION OR EVAPORATION OF QC  ******
         CND(I)=MAX(-QC(I), CND(I))
CC    ******   DEPOSITION OR SUBLIMATION OF QI    ******
         DEP(I)=MAX(-QI(I), DEP(I))
        PT(I)=PT(I)+AVCP*CND(I)+ASCP*DEP(I)
        QV(I)=QV(I)-CND(I)-DEP(I)
        QC(I)=QC(I)+CND(I)
        QI(I)=QI(I)+DEP(I)
  250   CONTINUE
C* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
C* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
        DO 275 I=2,ILES
         PSDEP(I)=0.0
         PSSUB(I)=0.0
         PGSUB(I)=0.0
          TAIR(I)=(PT(I)+TB0)*PI0
          IF(TAIR(I).LT.T0) THEN
             IF(QS(I).LT.cmin) QS(I)=0.0
             IF(QG(I).LT.cmin) QG(I)=0.0
             RTAIR(I)=1./(TAIR(I)-C76)
            QSI(I)=RP0*EXP(C218-C580*RTAIR(I))
            SSI(I)=(QV(I)+QB0)/QSI(I)-1.
              Y1(I)=R10AR/(TCA(I)*TAIR(I)**2)+1./(DWV(I)*QSI(I))
              Y2(I)=.78/ZS(I)**2+R101F*SCV(I)/ZS(I)**BSH5
           PSDEP(I)=R10T*SSI(I)*Y2(I)/Y1(I)
           PSSUB(I)=PSDEP(I)
            PSDEP(I)=MAX(PSDEP(I), 0.)
            PSSUB(I)=MAX(-QS(I), MIN(PSSUB(I), 0.))
              Y2(I)=.78/ZG(I)**2+R20BQ*SCV(I)/ZG(I)**BGH5
           PGSUB(I)=R20T*SSI(I)*Y2(I)/Y1(I)
            DM(I)=QV(I)+QB0-QSI(I)
            RSUB1(I)=CS580*QSI(I)*RTAIR(I)*RTAIR(I)
C     ********   DEPOSITION OR SUBLIMATION OF QS  **********************
              Y1(I)=DM(I)/(1.+RSUB1(I))
           PSDEP(I)=MIN(PSDEP(I),MAX(Y1(I),0.))
              Y2(I)=MIN(Y1(I),0.)
           PSSUB(I)=MAX(PSSUB(I),Y2(I))
C     ********   SUBLIMATION OF QG   ***********************************
              DD(I)=MAX((-Y2(I)-QS(I)), 0.)
           PGSUB(I)=MIN(DD(I), QG(I), MAX(PGSUB(I),0.))
              DLT1(I)=CVMGP(1.,0.,QC(I)+QI(I)-1.E-5)
            PSDEP(I)=DLT1(I)*PSDEP(I)
            PSSUB(I)=(1.-DLT1(I))*PSSUB(I)
            PGSUB(I)=(1.-DLT1(I))*PGSUB(I)
           PT(I)=PT(I)+ASCP*(PSDEP(I)+PSSUB(I)-PGSUB(I))
           QV(I)=QV(I)+PGSUB(I)-PSSUB(I)-PSDEP(I)
           QS(I)=QS(I)+PSDEP(I)+PSSUB(I)
           QG(I)=QG(I)-PGSUB(I)
          ENDIF
C* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
         ERN(I)=0.0
          IF(QR(I).GT.0.0) THEN
            TAIR(I)=(PT(I)+TB0)*PI0
            RTAIR(I)=1./(TAIR(I)-C358)
            QSW(I)=RP0*EXP(C172-C409*RTAIR(I))
            SSW(I)=(QV(I)+QB0)/QSW(I)-1.0
            DM(I)=QV(I)+QB0-QSW(I)
            RSUB1(I)=CV409*QSW(I)*RTAIR(I)*RTAIR(I)
             DD1(I)=MAX(-DM(I)/(1.+RSUB1(I)), 0.0)
CM             Y1(I)=R00*QR(I)
CM         ERN(I)=(((1.6+124.9*Y1(I)**.2046)*Y1(I)**.525)
CM   1          /(2.55E6/(P00*QSW(I))+5.4E5))*(-DM(I)/(R00*QSW(I)))*D2T
             Y1(I)=.78/ZR(I)**2+R23AF*SCV(I)/ZR(I)**BWH5
             Y2(I)=R23BR/(TCA(I)*TAIR(I)**2)+1./(DWV(I)*QSW(I))
           ERN(I)=R23T*SSW(I)*Y1(I)/Y2(I)
           ERN(I)=MIN(DD1(I),QR(I),MAX(ERN(I),0.))
           PT(I)=PT(I)-AVCP*ERN(I)
           QV(I)=QV(I)+ERN(I)
           QR(I)=QR(I)-ERN(I)
          ENDIF
C        IF (QV(I)+QB0 .LE. 0.) QV(I)=-QB0
        IF (QC(I) .LE. cmin) QC(I)=0.
        IF (QR(I) .LE. cmin) QR(I)=0.
        IF (QI(I) .LE. cmin) QI(I)=0.
        IF (QS(I) .LE. cmin) QS(I)=0.
        IF (QG(I) .LE. cmin) QG(I)=0.
C        Y1(I)=PT(I)-XX0(I)
C        Y2(I)=QV(I)-XX00(I)
C       DPT(I,K)=PT(I)+TMD(I,K)
C       DQV(I,K)=QV(I)+QMD(I,K)
C        TMD(I,K)=Y1(I)
C        QMD(I,K)=Y2(I)
        DPT(I,K)=PT(I)
        DQV(I,K)=QV(I)
        QCL(I,K)=QC(I)
        QRN(I,K)=QR(I)
        QCI(I,K)=QI(I)
        QCS(I,K)=QS(I)
        QCG(I,K)=QG(I)
  275   CONTINUE
CC*********************************************************************
       SCC=0.
       SEE=0.
       DO 110 I=2,ILES
         DD(I)=MAX(-CND(I), 0.)
         CND(I)=MAX(CND(I), 0.)
         DD1(I)=MAX(-DEP(I), 0.)
         DEP(I)=MAX(DEP(I), 0.)
         SCC=SCC+CND(I)
  110    SEE=SEE+DD(I)+ERN(I)
        SC(K)=SCC+SC(K)
        SE(K)=SEE+SE(K)
CC    ***   STATISTICS FOR CONVECTIVE AND ANVIL REGIMES   ***********
       IF (ID .EQ. 1) THEN
c         RDTS=1./D2T
        DO 280 I=2,ILES
c           CND(I)=CND(I)*RDTS
c           DD(I)=DD(I)*RDTS
c           PINT(I)=PINT(I)*RDTS
c           PIDW(I)=PIDW(I)*RDTS
c           PIMLT(I)=PIMLT(I)*RDTS
c           PIHOM(I)=PIHOM(I)*RDTS
c           PSMLT(I)=PSMLT(I)*RDTS
c           PGMLT(I)=PGMLT(I)*RDTS
c           PSDEP(I)=PSDEP(I)*RDTS
c           PSSUB(I)=PSSUB(I)*RDTS
c           PGSUB(I)=PGSUB(I)*RDTS
c           DD1(I)=DD1(I)*RDTS
c           DEP(I)=DEP(I)*RDTS
c           ERN(I)=ERN(I)*RDTS
         DDA(I)=RSW(I,K)
         DDB(I)=RLW(I,K)
         Y1(I)=QC(I)+QR(I)+QI(I)+QS(I)+QG(I)
         dm(i)=a0*(rho1(k)*ww1(i,k)+rho1(k+1)*ww1(i,k+1)+
     1              Y0(I)*(rho1(k)*wb(k)+rho1(k+1)*wb(k+1)))
  280    RQ(I)=.005*(RHO1(K)*(WW1(I,K)+wb(k))+
     1               RHO1(K+1)*(WW1(I,K+1)+wb(k+1)))/r00


        DO 1050 KC=1,7
c        KC=4

         DO 320 MT=1,4
         DO 320 I=2,ILES
          IBZ(I,MT)=0
  320      IF(ICS(I,mt).EQ.1) IBZ(I,MT)=1
        DO 315 I=2,ILES
            IBZ(I,1)=1
  315    CONTINUE
         DO 330 MT=1,4
         DO 330 I=2,ILES
           IF(KC.EQ.4) GO TO 330
           IF(KC.LE.3) GO TO 36
            IF (RQ(I).GT.RBY(KC)) IBZ(I,MT)=0
            GO TO 330
   36       IF (RQ(I).LT.RBY(KC)) IBZ(I,MT)=0
  330     CONTINUE
         DO 350 MT=1,4
          SWW=0.0
          SCC=0.0
          SEE=0.0
          A1=0.0
          A2=0.0
          A3=0.0
          A4=0.0
          A5=0.0
          A6=0.0
          A7=0.0
          A8=0.0
          A9=0.0
          A10=0.0
          A11=0.0
          A12=0.0
          A13=0.0
          A14=0.0
          A15=0.0
          A16=0.0
          A17=0.0
          A18=0.0
          A19=0.0
          A20=0.0
          A21=0.0
          A22=0.0
          A23=0.0
          A24=0.0
          A25=0.0
          A26=0.0
          A27=0.0
          A28=0.0
          A29=0.0
          A30=0.0
          A31=0.0
          A32=0.0
          A33=0.0
          A34=0.0
          A35=0.0
          A36=0.0
          A37=0.0
          A38=0.0
          A39=0.0
          A40=0.0
          A41=0.0
          A42=0.0
          A43=0.0
          A44=0.0
          A45=0
          A46=0.
          A47=0.
          A48=0.
          A49=0.
          A50=0.

         DO 30 I=2,ILES
          IF(IBZ(I,MT).EQ.1) THEN
           SWW=SWW+DM(I)
           SCC=SCC+CND(I)
           SEE=SEE+DD(I)
           A1=A1+PINT(I)
           A2=A2+(PIDW(I)+PIHOM(I))
           A3=A3+PIMLT(I)
           A4=A4+PSAUT(I)
           A5=A5+PSACI(I)
           A6=A6+PSACW(I)
           A7=A7+PRACI(I)
           A8=A8+PIACR(I)
           A9=A9+PRAUT(I)
           A10=A10+PRACW(I)
           A11=A11+PSFW(I)
           A12=A12+PSFI(I)
           A13=A13+(PGACS(I)+DGACS(I)+WGACS(I))
           A14=A14+DGACW(I)
           A15=A15+(DGACI(I)+WGACI(I))
           A16=A16+DGACR(I)
           A17=A17+(WGACR(I)+WGACI(I)+WGACS(I)+DGACW(I))
           A18=A18+PGAUT(I)
           A19=A19+PRACS(I)
           A20=A20+PSACR(I)
           A21=A21+PGFR(I)
           A22=A22+PSMLT(I)
           A23=A23+PGMLT(I)
           A24=A24+PSDEP(I)
           A25=A25+PSSUB(I)
           A26=A26+PGSUB(I)
           A27=A27+DD1(I)
           A28=A28+PRACI(I)*DLT3(I)
           A29=A29+PIACR(I)*DLT3(I)
           A30=A30+PSACR(I)*DLT2(I)
           A31=A31+DEP(I)
           A32=A32+(DGACR(I)+DGACI(I)+DGACS(I)+DGACW(I))
           A33=A33+QC(I)
           A34=A34+QR(I)
           A35=A35+QI(I)
           A36=A36+QS(I)
           A37=A37+QG(I)
           A38=A38+ERN(I)
           A39=A39+WGACR(I)
           A40=A40+QSACW(I)
          A41=A41+DDA(I)
          A42=A42+DDB(I)
            A43=A43+(QV(I)+QA1(K)-QB(K))
            A44=A44+(PT(I)+TA1(K)-TB(K))
            A45=A45+1.
            A46=A46+Y1(I)
            A47=a47+(PSACW(I)+PSFW(I)+DGACW(I)+PIACR(I)+DGACR(I)
     1               +WGACR(I)+PSACR(I)+PGFR(I)-PIMLT(I)+PIDW(I))
           A48=A48+(Y1(I)-QCL1(I,K)-QRN1(I,K)-QCI1(I,K)-QCS1(I,K)
     1            -QCG1(I,K))
           A49=A49+(QV(I)-DQV1(I,K))
           A50=A50+(PT(I)-DPT1(I,K))

          ENDIF
   30   CONTINUE
        SMF0(K,MT,KC)=SWW+SMF0(K,MT,KC)
        COC(K,MT,KC)=SCC+COC(K,MT,KC)
        COE(K,MT,KC)=SEE+COE(K,MT,KC)
        THOM(K,MT,KC)=THOM(K,MT,KC)+A1
        TDW(K,MT,KC)=TDW(K,MT,KC)+A2
        TMLT(K,MT,KC)=TMLT(K,MT,KC)+A3
        SAUT(K,MT,KC)=SAUT(K,MT,KC)+A4
        SACI(K,MT,KC)=SACI(K,MT,KC)+A5
        SACW(K,MT,KC)=SACW(K,MT,KC)+A6
        RACI(K,MT,KC)=RACI(K,MT,KC)+A7
        TACR(K,MT,KC)=TACR(K,MT,KC)+A8
        RAUT(K,MT,KC)=RAUT(K,MT,KC)+A9
        RACW(K,MT,KC)=RACW(K,MT,KC)+A10
        SFW(K,MT,KC)=SFW(K,MT,KC)+A11
        SFI(K,MT,KC)=SFI(K,MT,KC)+A12
        GACS(K,MT,KC)=GACS(K,MT,KC)+A13
        GACW(K,MT,KC)=GACW(K,MT,KC)+A14
        GACI(K,MT,KC)=GACI(K,MT,KC)+A15
        GACR(K,MT,KC)=GACR(K,MT,KC)+A16
        GWET(K,MT,KC)=GWET(K,MT,KC)+A17
        GAUT(K,MT,KC)=GAUT(K,MT,KC)+A18
        RACS(K,MT,KC)=RACS(K,MT,KC)+A19
        SACR(K,MT,KC)=SACR(K,MT,KC)+A20
        GFR(K,MT,KC)=GFR(K,MT,KC)+A21
        SMLT(K,MT,KC)=SMLT(K,MT,KC)+A22
        GMLT(K,MT,KC)=GMLT(K,MT,KC)+A23
        SDEP(K,MT,KC)=SDEP(K,MT,KC)+A24
        SSUB(K,MT,KC)=SSUB(K,MT,KC)+A25
        GSUB(K,MT,KC)=GSUB(K,MT,KC)+A26
        PERN(K,MT,KC)=PERN(K,MT,KC)+A27
        D3RI(K,MT,KC)=D3RI(K,MT,KC)+A28
        D3IR(K,MT,KC)=D3IR(K,MT,KC)+A29
        D2SR(K,MT,KC)=D2SR(K,MT,KC)+A30
        D2RS(K,MT,KC)=D2RS(K,MT,KC)+A31
        GDRY(K,MT,KC)=GDRY(K,MT,KC)+A32
        ERNS(K,MT,KC)=ERNS(K,MT,KC)+A38
        WGRS(K,MT,KC)=WGRS(K,MT,KC)+A39
        QSWS(K,MT,KC)=QSWS(K,MT,KC)+A40
        SRSW(K,MT,KC)=SRSW(K,MT,KC)+A41
        SRLW(K,MT,KC)=SRLW(K,MT,KC)+A42

        QC0(K,MT,KC)=A33
        QR0(K,MT,KC)=A34
        QI0(K,MT,KC)=A35
        QS0(K,MT,KC)=A36
        QG0(K,MT,KC)=A37
        SQC0(K,MT,KC)=SQC0(K,MT,KC)+QC0(K,MT,KC)
        SQR0(K,MT,KC)=SQR0(K,MT,KC)+QR0(K,MT,KC)
        SQI0(K,MT,KC)=SQI0(K,MT,KC)+QI0(K,MT,KC)
        SQS0(K,MT,KC)=SQS0(K,MT,KC)+QS0(K,MT,KC)
        SQG0(K,MT,KC)=SQG0(K,MT,KC)+QG0(K,MT,KC)

        qv0l(k,mt,kc)=a43
        tt0l(k,mt,kc)=a44
        sgpt(k,mt,kc)=a45
        tsqq(k,mt,kc)=a46
        stt0l(k,mt,kc)=stt0l(k,mt,kc)+tt0l(k,mt,kc)
        sqv0l(k,mt,kc)=sqv0l(k,mt,kc)+qv0l(k,mt,kc)
        ssgpt(k,mt,kc)=ssgpt(k,mt,kc)+sgpt(k,mt,kc)
        tsqq1(k,mt,kc)=tsqq1(k,mt,kc)+tsqq(k,mt,kc)

        q1tl(k,mt,kc)=q1tl(k,mt,kc)+a47
        sqqdtl(k,mt,kc)=sqqdtl(k,mt,kc)+a48
        sqdtl(k,mt,kc)=sqdtl(k,mt,kc)+a49
        stdtl(k,mt,kc)=stdtl(k,mt,kc)+a50

  350   CONTINUE
 1050  CONTINUE
cc    ********************************************
c      do 40 i=2,iles
c       cnd(i)=cnd(i)+dep(i)
c       ern(i)=ern(i)+dd(i)+dd1(i)
c       cnd(i)=cnd(i)*rft
c       ern(i)=ern(i)*rft
c       ww=rq(i)
c        iww=max(ww,-5.E0)
c        iww=min(ww,15.E0)
c       if(ww.ge.0.) scu1(k)=scu1(k)+cnd(i)
c       if(ww.lt.0.) sed1(k)=sed1(k)+ern(i)
c       if(iww.lt.1) go to 45
c        s9(iww,k)=s9(iww,k)+cnd(i)
c        s10(iww,k)=s10(iww,k)+ern(i)
c       go to 40
c  45   jww=iabs(iww)
c       if(jww.lt.1) go to 40
c        sn9(jww,k)=sn9(jww,k)+cnd(i)
c        sn10(jww,k)=sn10(jww,k)+ern(i)
c  40   continue

C--------------------------------------------------------------------
C    CONDENSATION:  cnd(ij)
C    EVAPORATION:   dd(ij)+ern(ij)
C    DEPOSITION:    dep(ij)+psdep(ij)+pint(ij)
C    SUBLIMATION:   dd1(ij)+pgsub(ij)-pssub(ij)
C    MELTING:       psmlt(ij)+pgmlt(ij)+pimlt(ij)
C    FREEZING:      pihom(ij)+pidw(ij)+psacw(ij)+psfw(ij)+dgacw(ij)+
c                   piacr(ij)+dgacr(ij)+wgacr(ij)+psacr(ij)+pgfr(ij)                
C    MASS FLUX:     DM(IJ)
C    CLOUD WATER:   QC(IJ)
C    RAIN:          
C    CLOUD ICE
C    SNOW
C    HAIL/GRAUPEL:
C----------------------------------------------------------------------
       DO 42 IJ=2,iles
        CND(IJ)=cnd(ij)*RFT
        ERN(IJ)=(ern(ij)+dd(ij))*RFT
        Y1(IJ)=(dep(ij)+psdep(ij)+pint(ij))*RFT
        Y2(IJ)=(dd1(ij)+pgsub(ij)-pssub(ij))*RFT
        Y3(IJ)=(psmlt(ij)+pgmlt(ij)+pimlt(ij))*RFT
        Y4(IJ)=(pihom(ij)+pidw(ij)+psacw(ij)+psfw(ij)+dgacw(ij)+
     1          piacr(ij)+dgacr(ij)+wgacr(ij)+psacr(ij)+pgfr(ij))*RFT
        Y5(IJ)=DM(IJ)
        A1=1.E6*R00*QR(IJ)                                                  
        A2=1.E6*R00*QS(IJ)                                                  
        A3=1.E6*R00*QG(IJ)                                                  
        A11=UCOR*(MAX(1.E-5,A1))**1.75
        A22=UCOS*(MAX(1.E-5,A2))**1.75
        A33=UCOG*(MAX(1.E-5,A3))**1.75
        ZDRY=MAX(1.,A11+A22+A33)
        DBZ(IJ)=10.*LOG10(ZDRY)
        IF (TAIR(IJ).GE.273.16) THEN
          A44=A11+UWET*(A22+A33)**.95
          ZWET=MAX(1.,A44)
          DBZ(IJ)=10.*LOG10(ZWET)
        ENDIF
        DBZ(IJ)=DBZ(IJ)
        QC(IJ)=QC(IJ)
        QR(IJ)=QR(IJ)
        QI(IJ)=QI(IJ)
        QS(IJ)=QS(IJ)
        QG(IJ)=QG(IJ)

   42  CONTINUE
        DO 44 IJ=it1,it2
        IF(RQ(IJ) .GE. 0.) SCU1(K)=SCU1(K)+CND(IJ)
        IF(RQ(IJ) .LT. 0.) SED1(K)=SED1(K)+ERN(IJ)

   44  CONTINUE
       DO 40 IJ=2,iles

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          IWW=0
          JWW=0
          IF (RQ(IJ).GT.-0.5) IWW=MIN(RQ(IJ)+0.5,15.)+1
          IF (RQ(IJ).LE.-0.5) JWW=MAX(RQ(IJ)-0.5,-5.)
          JWW=IABS(JWW)

        IF (ICS(IJ,2) .EQ. 1 .AND. IWW .GE. 1) THEN
          S9(IWW,K)=S9(IWW,K)+CND(IJ)
          S10(IWW,K)=S10(IWW,K)+ERN(IJ)
          S11(IWW,K)=S11(IWW,K)+Y1(IJ)
          S12(IWW,K)=S12(IWW,K)+Y2(IJ)
          S13(IWW,K)=S13(IWW,K)+Y3(IJ)
          S14(IWW,K)=S14(IWW,K)+Y4(IJ)
          S15(IWW,K)=S15(IWW,K)+Y5(IJ)
          S16(IWW,K)=S16(IWW,K)+QC(IJ)
          S17(IWW,K)=S17(IWW,K)+QR(IJ)
          S18(IWW,K)=S18(IWW,K)+QI(IJ)
          S19(IWW,K)=S19(IWW,K)+QS(IJ)
          S20(IWW,K)=S20(IWW,K)+QG(IJ)
          S21(IWW,K)=S21(IWW,K)+DBZ(IJ)
          SCNT(IWW,K)=SCNT(IWW,K)+1.

        ENDIF
        IF (ICS(IJ,2) .EQ. 1 .AND. JWW .GE. 1) THEN
          SN9(jww,K)=SN9(jww,K)+CND(IJ)
          SN10(jww,K)=SN10(jww,K)+ERN(IJ)
          SN11(jww,K)=SN11(jww,K)+Y1(IJ)
          SN12(jww,K)=SN12(jww,K)+Y2(IJ)
          SN13(jww,K)=SN13(jww,K)+Y3(IJ)
          SN14(jww,K)=SN14(jww,K)+Y4(IJ)
          SN15(jww,K)=SN15(jww,K)+Y5(IJ)
          SN16(jww,K)=SN16(jww,K)+QC(IJ)
          SN17(jww,K)=SN17(jww,K)+QR(IJ)
          SN18(jww,K)=SN18(jww,K)+QI(IJ)
          SN19(jww,K)=SN19(jww,K)+QS(IJ)
          SN20(jww,K)=SN20(jww,K)+QG(IJ)
          SN21(jww,K)=SN21(jww,K)+DBZ(IJ)
          SNCNT(jww,K)=SNCNT(jww,K)+1.

        ENDIF
        IF (ICS(IJ,3) .EQ. 1 .AND. IWW .GE. 1) THEN
          SS9(IWW,K)=SS9(IWW,K)+CND(IJ)
          SS10(IWW,K)=SS10(IWW,K)+ERN(IJ)
          SS11(IWW,K)=SS11(IWW,K)+Y1(IJ)
          SS12(IWW,K)=SS12(IWW,K)+Y2(IJ)
          SS13(IWW,K)=SS13(IWW,K)+Y3(IJ)
          SS14(IWW,K)=SS14(IWW,K)+Y4(IJ)
          SS15(IWW,K)=SS15(IWW,K)+Y5(IJ)
          SS16(IWW,K)=SS16(IWW,K)+QC(IJ)
          SS17(IWW,K)=SS17(IWW,K)+QR(IJ)
          SS18(IWW,K)=SS18(IWW,K)+QI(IJ)
          SS19(IWW,K)=SS19(IWW,K)+QS(IJ)
          SS20(IWW,K)=SS20(IWW,K)+QG(IJ)
          SS21(IWW,K)=SS21(IWW,K)+DBZ(IJ)
          SSCNT(IWW,K)=SSCNT(IWW,K)+1.

        ENDIF
        IF(ICS(IJ,3) .EQ. 1 .AND. jww.GE.1) THEN
          SSN9(jww,K)=SSN9(jww,K)+CND(IJ)
          SSN10(jww,K)=SSN10(jww,K)+ERN(IJ)
          SSN11(jww,K)=SSN11(jww,K)+Y1(IJ)
          SSN12(jww,K)=SSN12(jww,K)+Y2(IJ)
          SSN13(jww,K)=SSN13(jww,K)+Y3(IJ)
          SSN14(jww,K)=SSN14(jww,K)+Y4(IJ)
          SSN15(jww,K)=SSN15(jww,K)+Y5(IJ)
          SSN16(jww,K)=SSN16(jww,K)+QC(IJ)
          SSN17(jww,K)=SSN17(jww,K)+QR(IJ)
          SSN18(jww,K)=SSN18(jww,K)+QI(IJ)
          SSN19(jww,K)=SSN19(jww,K)+QS(IJ)
          SSN20(jww,K)=SSN20(jww,K)+QG(IJ)
          SSN21(jww,K)=SSN21(jww,K)+DBZ(IJ)
          SSNCNT(jww,K)=SSNCNT(jww,K)+1.

        ENDIF
   40  CONTINUE

       ENDIF
 1000 CONTINUE

c     ****************************************************************
      if (id .eq. 1) then
       lc=0
       ls=0
       ln=0
       do 390 i=2,iles
        lc=lc+ics(i,2)
        ls=ls+ics(i,3)
        ln=ln+ics(i,4)
  390  continue
        if (lc .eq. 0) lc=1000000
        if (ls .eq. 0) ls=1000000
        if (ln .eq. 0) ln=1000000
       do 400 mt=1,4
        a1=ril2
         if (mt .eq. 2) a1=1./float(lc)
         if (mt .eq. 3) a1=1./float(ls)
         if (mt .eq. 4) a1=1./float(ln)
       do 410 k=2,kles
        y1(k)=0.0
        y2(k)=0.0
       do 410 i=2,iles
        if(ics(i,mt).eq.1) then
         y1(k)=y1(k)+dpt(i,k)
         y2(k)=y2(k)+dqv(i,k)
        endif
  410  continue
       do 430 k=2,kles
         tbt(k,mt)=y1(k)*a1
  430    qbt(k,mt)=y2(k)*a1
  400 continue
      endif
c     ****************************************************************
c     ****************************************************************
      do 460 k=2,kles
       y1(k)=0.
       y2(k)=0.
       sq(k)=0.
       sqc(k)=0.
       sqr(k)=0.
       sqi(k)=0.
       sqs(k)=0.
  460  sqg(k)=0.
      do 465 i=2,iles
      do 465 k=2,kles
       y1(k)=y1(k)+dpt(i,k)
       y2(k)=y2(k)+dqv(i,k)
       sq(k)=sq(k)+qcl(i,k)+qrn(i,k)+qci(i,k)+qcg(i,k)+qcs(i,k)
       sqc(k)=sqc(k)+qcl(i,k)
       sqr(k)=sqr(k)+qrn(i,k)
       sqi(k)=sqi(k)+qci(i,k)
       sqs(k)=sqs(k)+qcs(i,k)
       sqg(k)=sqg(k)+qcg(i,k)
  465 continue
      do 470 k=2,kles
       y1(k)=y1(k)*ril2
       y2(k)=y2(k)*ril2
       sq(k)=sq(k)*ril2
       sqc(k)=sqc(k)*ril2
       sqr(k)=sqr(k)*ril2
       sqi(k)=sqi(k)*ril2
       sqs(k)=sqs(k)*ril2
       sqg(k)=sqg(k)*ril2
       stqc(k)=stqc(k)+sqc(k)
       stqr(k)=stqr(k)+sqr(k)
       stqi(k)=stqi(k)+sqi(k)
       stqs(k)=stqs(k)+sqs(k)
       stqg(k)=stqg(k)+sqg(k)
  470 continue
      do 500 k=2,kles
       do 510 i=2,iles
        dpt(i,k)=dpt(i,k)-y1(k)
        dqv(i,k)=dqv(i,k)-y2(k)
  510  continue
       if(n.gt.1) go to 515
         IF (IJKADV .EQ. 1) THEN
           Y3(K)=TA(K)
           Y4(K)=QA(K)
          TA(K)=TA(K)+Y1(K)
          QA(K)=QA(K)+Y2(K)
          TA1(K)=y3(K)
          QA1(K)=y4(K)
         ELSE
          TA(K)=TA1(K)+Y1(K)
          QA(K)=QA1(K)+Y2(K)
         ENDIF
       go to 516
  515  do 517 i=2,iles
        IF (IJKADV .EQ. 0) THEN
         DPT1(I,K)=DPT1(I,K)+EPS*DPT(I,K)
         DQV1(I,K)=DQV1(I,K)+EPS*DQV(I,K)
         QCL1(I,K)=QCL1(I,K)+EPS*QCL(I,K)
         QRN1(I,K)=QRN1(I,K)+EPS*QRN(I,K)
         QCI1(I,K)=QCI1(I,K)+EPS*QCI(I,K)
         QCS1(I,K)=QCS1(I,K)+EPS*QCS(I,K)
         QCG1(I,K)=QCG1(I,K)+EPS*QCG(I,K)
c        ELSE
c          DPT1(I,K)=DPT(I,K)
c          DQV1(I,K)=DQV(I,K)
c          QCL1(I,K)=QCL(I,K)
c          QRN1(I,K)=QRN(I,K)
c          QCI1(I,K)=QCI(I,K)
c          QCS1(I,K)=QCS(I,K)
c          QCG1(I,K)=QCG(I,K)
        ENDIF
  517  continue
         IF (IJKADV .EQ. 1) THEN
           Y3(K)=TA(K)
           Y4(K)=QA(K)
          TA(K)=TA(K)+Y1(K)
          QA(K)=QA(K)+Y2(K)
          TA1(K)=y3(K)
          QA1(K)=y4(K)
         ELSE
           A0=TA(K)+EPS*(-2.*TA(K)+TA1(K))
          TA(K)=TA1(K)+Y1(K)
          TA1(K)=A0+EPS*TA(K)
           A0=QA(K)+EPS*(-2.*QA(K)+QA1(K))
          QA(K)=QA1(K)+Y2(K)
          QA1(K)=A0+EPS*QA(K)
         ENDIF
  516  st(k)=ta(k)-tb(k)
       sv(k)=qa(k)-qb(k)
  500 continue
c     ****   fix the boundary conditions
      call boundy (dpt,dpt1)
      call boundy (dqv,dqv1)
      call boundy (qcl,qcl1)
      call boundy (qrn,qrn1)
      call boundy (qci,qci1)
      call boundy (qcs,qcs1)
      call boundy (qcg,qcg1)
c     ***********************************
      sq(1)=sq(2)
      ta(1)=ta(2)
      qa(1)=qa(2)
      ta1(1)=ta1(2)
      qa1(1)=qa1(2)
      sq(kmax)=sq(kles)
      ta(kmax)=ta(kles)*2.-ta(kl2)
      qa(kmax)=qa(kles)*2.-qa(kl2)
      ta1(kmax)=ta1(kles)*2.-ta1(kl2)
      qa1(kmax)=qa1(kles)*2.-qa1(kl2)

      do 600 k=3,kles
       IF (IJKADV .EQ. 1) THEN
        FD(K)=(TA1(K)-TA1(K-1))*RHO1(K)*.5
        FE(K)=(QA1(K)-QA1(K-1))*RHO1(K)*.5
       ELSE
        FD(K)=(TA(K)-TA(K-1))*RHO1(K)*.5
        FE(K)=(QA(K)-QA(K-1))*RHO1(K)*.5
       ENDIF
  600 continue 
      fd(kmax)=fd(kles)
      fe(kmax)=fe(kles)
      d2t=d22t
      print*, 'right here'
      stop
      RETURN
      END
