c              IF USING CRAY, COMPILE WITH -Ez FLAG 
C              uncomment following line if using sgi, add comment if using cray!

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                                                                      c
      PROGRAM COARE_GATE_3D
c                                                                      c
c    *****   7-21-99                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ******************************************************************
C     **   ANELASTIC VERSION CYCLIC LATERAL BOUNDARY CONDITIONS       **
C     **   filename=coare-gate3d.f           (STARTS 4-4-97)          **
C     **   TOGA COARE FLUX ALGORITHM   (STARTS 6-6-95)                **
C     **   LONGWAVE AND SHORTWAVE RADIATION (STARTS 6-6-95)           **
C     **   POSITIVE DEFINITE ADVECTION (STARTS 2-27-94)               **
C     **   3-D RUN CRAY/YMP VERSIONN (STARTS 11/09/90)                **
C     **   3-D RUN WITH HALF-PRECISION (STARTS 6/24/84)               **
C     **   ICE-PHASE (STARTS 10/24/84, 4/12/86)                       **
C     **   SECOND ORDER FLUX FORM (STARTS 3/9/78)                     **
C     **      BOB"S DIFFUSION COEFFICIENT                             **
C     **   THREE-D AND VECTORIZED FORM FOR CRAY-1 (STARTS 9/12/80)    **
C     **      STRETCHED IN Z-DIRECTION                                **
C     **      STRETCHED IN X-DIRECTION        (STARTS 3/21/89)        **
C     **      OPEN LATERAL B.C. IN X-AXIS     (STARTS 3/21/89)        **
C     **      PERIODIC L. B. C. IN Y-AXIS     (STARTS 3/21/89)        **
C     **      COMPRESSIBLE SYSTEM             (STARTS 3/21/89)        **
C     **   FOR CYBER 205 VECTORIZED PROCESSOR (STARTS 12/18/82)       **
C     ******************************************************************
cccshie 10/23/01, NT=2880 for 10-day rain with imx=300(sec)
cccshie 11/16/01, imx=900(sec)=15min (control print max values, nothing to do with
ccc               rain in 3d ("imx" control rain output in 2d)

c     PARAMETER (NX=258,NY=258,NZ=43,NM=NX,NT 5760,ITT=244,lay=88) ! dan's original

      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880,ITT=244,lay=88) ! shie 10/23/01

      PARAMETER (NXI=NX-2,NYJ=NY-2,NT2=2*NT,NM12=12*NM,itt3=3*itt)

       common /setbound/ icheck

      common/bch1/ rby(7)
      COMMON/SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
      COMMON/SURFACE/ SUN_4(NX,NY,4)
      COMMON/BPBL/ UHT(NZ),WHT(NZ),TGBAT0
      COMMON/BNDY2/ UB0(NZ),VB0(NZ)
      COMMON/BGWY/ UBAR,VBAR
      COMMON/ISFCRI/ IRICONT

      common/sfcri/ ri180(nx,ny),riold180(nx,ny)
      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      common/dinradn/ qairsfc(nx,ny),pairsfc(nx,ny),thairsf(nx,ny)
      common/dinrad1/ dz1half
c
      COMMON/SLTM/ RLAT,RMONTH,RIDAY,HRL,SO0,COSZ,ICOSZ,TERMAN,RSFC
      common/srflx/ sfir(lay),sfne(lay),qswm(lay),qirm(lay)
      common/clear/ wave(nz),tave(nz),ts000,qs000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      common/msui/ isounding,isui,kfine
      common/option/ lipps,ijkadv,istatmin,iwater,itoga,imlifting,lin,
     1   irf,iadvh,irfg,ismg,id
      common/option1/ icomp,nsml
      common/option2/ isngtao,iwbar,iuvbar,isfc,ice,ice2,iradave,idelz
c
      common/timestat/ ndt_stat,idq

      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1   NRAN,KT1,KT2
      COMMON/SYMB/ BLK,DOT,COM,EXM,PICE,AST,AGO,ABT
      COMMON/RTIME/ ISEC1
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
c
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSKT2(NZ),BSKM(NZ),BSKM4(NZ),
     1   BSIT(NX,NY),BSIT2(NX,NY),BSIM(NX,NY),BSIM4(NX,NY),
     2   BSJT(NX,NY),BSJT2(NX,NY),BSJM(NX,NY),BSJM4(NX,NY)
      COMMON/B4Z/ BSITZ(NZ),BSIT2Z(NZ),BSIMZ(NZ),BSIM4Z(NZ),BSJTZ(NZ),
     1   BSJT2Z(NZ),BSJMZ(NZ),BSJM4Z(NZ)
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
c
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B2S/ QCS1(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
      COMMON/B4A/ AK1(NX,NY,NZ)
      COMMON/SLWAVE/ RSW(NX,NY,NZ),RLW(NX,NY,NZ)
      COMMON/BSAT/ XXW(NX,NY,NZ)
      COMMON/BSAT1/ AAA(NX,NY,NZ)
      COMMON/BADV/ DFC(NX,NY,NZ)
      COMMON/BADV1/ DF0(NX,NY,NZ)
      COMMON/B3CR/ QCL2(NX,NY,NZ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
      COMMON/BBXX/ XXX(NX,NY,NZ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      common/b66b/ s_dep(nz),s_sub(nz),s_qrs(nz),s_qrl(nz),s_mel(nz),
     1   s_frz(nz)
      COMMON/B7/ SQAAQ(NZ),SQTT(NZ),SQVV(NZ),SQAAT(NZ),SQAK(NZ)
      COMMON/BRH1/ SRRO(NZ),QRRO(NZ),SQC(NZ),SQR(NZ),SQI(NZ),SQS(NZ),
     1   SQG(NZ),STQC(NZ),STQR(NZ),STQI(NZ),STQS(NZ),STQG(NZ),
     2   TQC(NT),TQR(NT),TQI(NT),TQS(NT),TQG(NT)
      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      COMMON/BW1/ PCLTOP(NXI,NYJ),PCLBOT(NXI,NYJ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/B8/ SMF(NZ),SMU(NZ),SMD(NZ),STF(NZ),STU(NZ),STD(NZ),
     1   SQF(NZ),SQU(NZ),SQD(NZ),SDT(NZ),SDQ(NZ),CLD(NZ),STV(NZ)
      COMMON/B9/ SUT1(NZ),SUC1(NZ),SUN1(NZ),SUU1(NZ),SUD1(NZ),AUB(NZ),
     1   SSU(NZ),SUT2(NZ),SUC2(NZ),SUN2(NZ),SUU2(NZ),SUD2(NZ),AUC1(NZ),
     2   AUN1(NZ),AUU1(NZ),AUD1(NZ),AUC2(NZ),AUN2(NZ),AUU2(NZ),AUD2(NZ),
     3   SVT1(NZ),SVC1(NZ),SVN1(NZ),SVU1(NZ),SVD1(NZ),AVB(NZ),
     4   SSV(NZ),SVT2(NZ),SVC2(NZ),SVN2(NZ),SVU2(NZ),SVD2(NZ),AVC1(NZ),
     5   AVN1(NZ),AVU1(NZ),AVD1(NZ),AVC2(NZ),AVN2(NZ),AVU2(NZ),AVD2(NZ),
     6   STF1(NZ),STU1(NZ),STD1(NZ),SQF1(NZ),SQU1(NZ),SQD1(NZ),
     7   STF2(NZ),STU2(NZ),STD2(NZ),SQF2(NZ),SQU2(NZ),SQD2(NZ),CLDU(NZ),
     8   CLDD(NZ),CLDN(NZ),SMN(NZ),STN1(NZ),STN2(NZ),STN(NZ),SQN1(NZ),
     9   SQN2(NZ),SQN(NZ)
      COMMON/BSTS/ THOM(NZ,4,7),TDW(NZ,4,7),TMLT(NZ,4,7),SAUT(NZ,4,7),
     1 SACI(NZ,4,7),SACW(NZ,4,7),RACI(NZ,4,7),TACR(NZ,4,7),RAUT(NZ,4,7),
     2 RACW(NZ,4,7),SFW(NZ,4,7),SFI(NZ,4,7),GACS(NZ,4,7),GACW(NZ,4,7),
     3 GACI(NZ,4,7),GACR(NZ,4,7),GWET(NZ,4,7),GAUT(NZ,4,7),RACS(NZ,4,7),
     4 SACR(NZ,4,7),GFR(NZ,4,7),SMLT(NZ,4,7),GMLT(NZ,4,7),SDEP(NZ,4,7),
     5 SSUB(NZ,4,7),GSUB(NZ,4,7),PERN(NZ,4,7),D3RI(NZ,4,7),D3IR(NZ,4,7),
     6 D2SR(NZ,4,7),D2RS(NZ,4,7),GDRY(NZ,4,7),COC(NZ,4,7),COE(NZ,4,7),
     7 SMF0(NZ,4,7),QC0(NZ,4,7),QR0(NZ,4,7),QI0(NZ,4,7),QS0(NZ,4,7),
     8 QG0(NZ,4,7),SQC0(NZ,4,7),SQR0(NZ,4,7),SQI0(NZ,4,7),SQS0(NZ,4,7),
     9 SQG0(NZ,4,7),ERNS(NZ,4,7),WGRS(NZ,4,7),QSWS(NZ,4,7),TB0(NZ,4),
     1 QB0(NZ,4)
      COMMON/BSTS1/ TUT1(NZ,4,7),TUT2(NZ,4,7),TVT1(NZ,4,7),TVT2(NZ,4,7),
     1 TSTF(NZ,4,7),TSTF1(NZ,4,7),TSTF2(NZ,4,7),TSQF(NZ,4,7),
     2 QQQ(NZ,4,7),TSQF1(NZ,4,7),TSQF2(NZ,4,7),TSQQ(NZ,4,7),
     3 TSQQ1(NZ,4,7)
      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNTL(NZ,4,7),
     1 SNQH(NZ,4,7),SNQV(NZ,4,7),SNQD(NZ,4,7),SNQL(NZ,4,7),
     2 SNTL1(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7),
     3 SNQL1(NZ,4,7)
      COMMON/BSTS3/ QV0(NZ,4,7),TT0(NZ,4,7),SQV0(NZ,4,7),STT0(NZ,4,7),
     1 SGPT(NZ,4,7),SSGPT(NZ,4,7),SNQHD(NZ,4,7),SNQVD(NZ,4,7),
     2 Q1Tt(NZ,4,7),SNHDH(NZ,4,7),SQHDT(NZ,4,7),SQVDT(NZ,4,7)
      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7),SQTDT(NZ,4,7),SQHL(NZ,4,7)
      COMMON/BSTS20/ OTHERT_ADD(NZ,4,7),OTHERQ_ADD(NZ,4,7)

      common/bsts40/ fcld(nz,4,7)

      COMMON/BTT/ S1(16,NZ,14),SN1(5,NZ,14)
      COMMON/BCS/ S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     1 S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     2 S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     3 SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     4 SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     5 SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCSS/ SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     1 SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     2 SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     3 SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     4 SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     5 SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/BBA/ XYP(NX,NY),XY(NX,NY),XYM(NX,NY),XYPM(NX,NY),Y1(NM),
     1   Y2(NM),Y3(NM),Y4(NM),Y5(NM12)
      COMMON/BI/ IT(NX,NY),ICS(NX,NY,4)
      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NX,NY),ICS5(NX,NY,4),
     1  IBZ(NX,NY,4)
      COMMON/RSTAT/ CSTTT(NX,NY),CSTT(NX,NY)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
      common/q1q2z/ q1z_6h(nz,itt),q1zt(nz),q2z_6h(nz,itt),q2zt(nz)
      common/blsn/ ths(itt),qs(itt),ts(itt),
     1   pss(itt)

      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/gbs/ tlsw(nz),qlsw(nz),ttlsw(nz),tqlsw(nz),smfff(nz),
     1   smffff(nz),tcof(nz),textra(nz),qextra(nz)
      common/bbb6/ wb1(nz),wb2(nz),wb3(nz)
      common/bb6/ tls(nz),qls(nz),tls1(nz),tls2(nz),qls1(nz),qls2(nz),
     1   tls3(nz),tls4(nz),qls3(nz),qls4(nz),sft(nz),sfq(nz)
      common/sfcten/ tsdt,qsdt,thsdt,psdt
      common/TMP/ y1d(nx,ny,nz),y2d(nx,ny,nz),y3d(nx,ny,nz),
     1            y4d(nx,ny,nz),y03d(nx,ny,1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DIMENSION FV(NZ),FV1(NZ),QAA1(NZ),y6(nm)
      COMMON/BSFC/ tsfc_1(nx,ny), qsfc_1(nx,ny)
      common/bls3/ factor_nuding ! control nudging
      common/bls4a/ ubi(nz),vbi(nz),ub_2h(nz,itt3),vb_2h(nz,itt3)
      common/tvertical/ denave

cccshie 11/16/01 for 3d!
      common/iceopt/ ice913,ilif
      common/iice/ new_ice_sat

      integer ised(12)
      DATA RBY/2.,1.,0.,0.,0.,-.5,-1./
      data ised/54368,54368,33214,43598,21368,11246,78246,62359,31246,
     1   29548,83124,98272/
      character*4 wp,up,vp,pt,sh,cd,rn,pp,rs,rg,df,qi,shd
      character*1 blk,dot,com,exm,pice,ast,ago,abt
c     data wp,up,vp,pt,sh,cd,rn,pp,rs,rg,df,qi/'W   ','U   '

cccshie 11/16/01 add lable shd='DFO '
      data wp,up,vp,pt,sh,cd,rn,pp,rs,rg,df,qi,shd/'W   ','U   '
     1     ,'V   ','DPT ','DQV ','CLD ','RN  ','PI  '
     2     ,'QS  ','QG  ','AK  ','QI  ','DFO '/
      data blk,dot,com,exm,pice,ast,ago,abt/' ','.',',','R','I','*'
     1     ,'G','O'/
      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ICASE=1:  Sui et al. 94   iwbar=1
C     ICASE=0:  TOGA/GATE       iwbar=0
C     ICE = 0       WARM RAIN ONLY
C     ICE2=1        QI AND QS ONLY
C     ICERH=0, ICELIN=1:  FOR LIN-TYPE 3 CLASSES ICE SCHEME
C     ICERH=1, ICELIN=0:  FOR RH-TYPE 3 CLASSES ICE SCHEME
C     IJKADV=0       HIGHER ORDER (2ND OR 4TH ORDER) ADVECTION
C     IJKADV=1       LOWERER ORDERBUT POSITIVE DEFINITE ADVECTION
C     ILIF=0  (MODEL INTERCOMPARISON WORKSHOP)
C     IMLIFTING=1, WITH LARGE-SCALE LIFTING (SOONG & OGURA, TAO & SOONG)
C     IMLIFTING=0, NO LARGE-SCALE LIFTING (TRADITIONAL MODEL)
C     INCAR = 0 for Goddard data set - different unit than NCAR
C     IRAD = 0, NO RADIATIVE PROCESSES
C     IRADAVE=1:  Include cloud effects but apply Qr and sfc fluxes uniformly over d omain
C     ISFC=1:  INCLUDE SURFACE PROCESSES
C     ISFC=0:  NO SURFACE FLUXES
C     ISFC HAS TO BE 1 IN ORDER TO HAVE ITOGA=1 OR 0
C     ISNGTAO=0      WBAR * D (THEA)/DZ IS CONSTANT OR SPECIFIED
C     ISNGTAO=1      WBAR IS THE MAIN FORCING & LAPSE RATE VARIES SEE SOONG AND TAO 1980
C     ISOUNDING = 1 USE SUI ET AL'S SOUNDING
C     ISOUNDING = 0 USE GRABOWSKI ET AL'S SOUNDING
c     isui=1  following sui et al (1994)  WBAR * [D(THEA)/DZ]T
c     isui=0  following grabowski et al (1996) [WBAR * D(THEA)/DZ]CONST
C     ITEST=0  THERE BE NO STATISTICS
c     itest_model=1: use 913's ocean sfc t and qv
C     ITOGA=0:      BLACKADAR PBL
C     ITOGA=1:      TOGACOARE SURFACE FLUX ROURINE
C     ITOGA=2       AERODYNAMIC FORMULA
C     IUVBAR = 1    U/VBAR ARE AS IMPOSED-SPECIFIED SEE SOONG AND TAO 1980
C     IUVBAR = 0    U/VBAR VARY (NO CONTROL) SEE TAO AND SOONG 1986
C     IWBAR = 1     WBAR IS CONSTANT WITH TIME - SEE SUI ET AL 1994
C     IWBAR = 0     WBAR VARIES AS IMPOSED TAO AND SOONG  SEE SOONG AND TAO 1980
C     KFINE=0,  DZ(1)=200-250 M; SUI ET AL (1994)  COARSE VERTICAL RESOLUTION
C     KFINE=1 GRABOWSKI ET AL (1996) FINE VERTICAL RESOLUTION
C               USE FINE VERTICAL RESOLUTION AT FIRST VERTICAL GRID PT
C               THIS IS NECESSARY FOR RUN USING SURFACE FLUX ALGORITHM
C     LIN=1, HAIL CASE - MIDLATITUDEs
C     LIN=0, GRAUPEL CASE - TROPICs 
C     LIPPS = 1      FOLLOWED LIPPS & HEMLER'S MODIFICATION
C     LIPPS = 0      ORIGINAL ANELASTIC SYSTEM
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     it starts -6 h and all temporal forcing are zero -> t,q,u,v and w
C                    are constant
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ifmicro=15
      ifdyn=16
      ifstat=17
      ifrain=18
      open(ifmicro,file='data_micro',form='unformatted')
      open(ifdyn,file='data_dyn',form='unformatted')
      open(ifstat,file='data_stat',form='unformatted')
      open(ifrain,file='data_rain',form='unformatted')

      cpu=0.

C - Begin CPU timing of simulation

c     ibudsec=300
c     rbud = 1.

      call timing(time1)

C - IRS = 1 for restart.  Allocate memory and initialize all variables to 0
cccshie
      irs = 0
      call rinit (irs)

C  Time Flags
c
cccshie
      dt = 12.
      d2t = dt
cccshie
ccc   time=8.*24.*60.  ! in minutes, final goal: 8 days for may18-may26
ccc   time=9.*24.*60.  ! in minutes, final goal: 9 days for june2-june11

c     time =  0.*24.*60. + 1.*60.  ! 10/23/01 shie, 1hr trial run
      time =  0.*24.*60. + 3.*60.  ! 10/29/01, starting real run, every 3-hr

c     time =  0.*24.*60. + 6.*60.
c      time = 0.*24.*60.+ 12.*60.
c      time = 0.*24.*60.+ 18.*60.
c      time = 1.*24.*60.+ 0.*60.
      sec = dt
      aminut = sec / 60.
      iend = time * 60. + 1.
      DTS = 2.
      EPS = .1
      NSML = dt/dts

C  Microphysics Flags

      ice = 1
      ice2 = 0
      icelin = 0
      icerh = 1
      lin = 0
      iwater = 0
        if (ice .eq. 1) iwater = 0
        if (iwater .eq. 1) ice = 0
        if (icelin .eq. 1) icerh = 0
        if (icerh .eq. 1) icelin = 0


C  Grid Setup Flags

      icase = 0
      idelz = 0

      kdist = 2
      kfine = 0
      if (kfine .eq. 0) kdist=2

c      kdist=4
c       if(.01*z2(3).ge.300.) kdist=3
      print*,'kdist is',kdist

      fcor = 0.
      rlat = 0.
      fcorg = fcor / 980.
      icomp = 0

C  Advection Scheme Flags

      iadvh = 4
      iadvh1 = 4
      ijkadv = 1
      ismg = 0

C  Flags To Change How Forcing is Applied

cccshie
c     imlifting = 1
      imlifting = 0 ! 10/29/01, following 2d run

      isngtao = 0
      isounding = 0
      istrt = 21600
cccshie
c     lipps = 1
      lipps = 0 ! 10/29/01, following 2d run
      isui = 0
      iuvbar = 1
      iwbar = 0

      ir36 = 21600
      r36 = 1. / float(ir36)

      incar = 1
cccshie 11/16/01 for 3d so far, it is only true for "ilif = 0", more modifications
ccc      are needed for "ilif = 1" to work (refernce: "scsmex 2d code")
      ilif = 0

cccshie by tao 5/3/01, shie 11/16/01 3d
      new_ice_sat = 1 !  call tao et al (1989) saturation technique twice
c     new_ice_sat=0
c


      irf = 1
      irfg = 0

      itest = 1
      itest_model = 0
     
C  Surface Flags

      isfc = 1
cccshie
      itoga = 2 ! 11/01/01, since "itoga = 1" did not work so far.
c     itoga = 1 ! 10/29/01, following 2d run
        if (isfc .eq. 0) itoga = 10
CTAO
c      isfcave = 1:  Use mean U and V for flux calculation
       isfcave=0
CTAO
      idpbl = 0
      pbldt = 180.
      ipbldt = pbldt
      iricont = 0
      iricont1 = 0

C  Radiation Flags

cccshie
      hrl = 8.0
cc    riday = 15. ! dan's original

      riday=18  ! 10/23/01 long term integ, 051800-052600 ( 0+18=18)
c     riday=33  !          long term integ, 060200-061100 (31+ 2=33)

      rmonth = 5.

      irad = 1
      iradave = 0
      iradd = 1
      ircl =  0
c     IRADDT = 180
      IRADDT = 600  ! 10/29/01, to speed up the model running
        if (irad .eq. 0) iradd = 0

C - Statistics and Printout Flags

cccshie
      INP = 3600
      IMP = 3600
c     IMX = 300  ! 5 minutes
      IMX = 900  ! cccshie 11/16/01 changed to 15 minutes, to print maximum value

      ISTAT  = 1800 ! shie 10/23/01
      IMICRO = 1800 ! shie 10/23/01
      IDYN   = 1800 ! shie 10/23/01
c     ISTAT = 900   ! dan's original
c     IMICRO = 900  ! dan's original
c     IDYN = 900    ! dan's original

      IWTP = 0
c
      ndt_stat=180
      ISTATMIN = ndt_stat
c
      CE=.7
      CK=.2
C
      ID=0
      M=0
      N=1
C     **************
      UBAR=0.
      VBAR=0.

      CALL BASE(IRS)
      IF (IWATER .EQ. 1) THEN
      ELSE
         IF (ICERH.EQ.1) CALL CONSATRH
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL SLVPI (0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call zang (cosz,rlat,rmonth,riday,hrl)
      ircl=0
c      iflag=0

      IF (IRAD .EQ. 1) CALL RADRAT (0)

c      iflag=1

      CCON=1.E-5
      NGTP=1
      IRMS=2
      IRMF=iles
      JRMS=2
      JRMF=jles
C     ******************************************************************
      WRITE(6,155)
      WRITE(6,154)'IL2=',IL2
      WRITE(6,154)'JL2=',JL2
      WRITE(6,154)'KL2=',KL2
      WRITE(6,1544)'DX=',DX
      WRITE(6,1544)'DY=',DY
      WRITE(6,1544)'DZ=',DZ
      WRITE(6,1544)'TIME=',TIME
      WRITE(6,1544)'DT=',DT
      WRITE(6,1544)'PSFC=',PSFC
      WRITE(6,1544)'RD1=',RD1
      WRITE(6,1544)'RD2=',RD2
      WRITE(6,1544)'BOUND=',BOUND
      WRITE(6,1544)'CK=',CK
      WRITE(6,1544)'CE=',CE
      WRITE(6,1544)'RM=',RM0
      WRITE(6,1544)'EPS=',EPS
      WRITE(6,154)'INP=',INP
      WRITE(6,154)'IMP=',IMP
      WRITE(6,154)'ISED=',ISED(2)
      WRITE(6,154)'NGTP=',NGTP
      WRITE(6,1544)'COTO=',CCON
      WRITE(6,1544)'RM1=',RM1

  154 FORMAT (A,I12)
 1544 FORMAT (A,F13.2)

C     ******************************************************************
       kt1=2
       kt2=kles
      DO 10 K=1,KMAX
       SFT(K)=.0
       SFQ(K)=.0
       QAA1(K)=0.
       FV(K)=SQRT(RHO(2)*RRHO(K))
       FV1(K)=SQRT(RHO1(2)*RRHO1(K))
       q1t(K)=0.
       q2t(K)=0.
       WBT(K)=0.
   10 CONTINUE

c$doacross local(j,i)
      DO J=1,JMAX
      DO I=1,IMAX
       SUW(I,J)=0.
       SVW(I,J)=0.
       SWT(I,J)=0.
       SWQ(I,J)=0.
       ri180(i,j)=0.
       riold180(i,j)=0.
       y0(i,j)=1.
       y03d(i,j,1)=1.
       CSTT(I,J)=0.
       CSTTT(I,J)=0.
       RX(I,J)=0.
      ENDDO
      ENDDO

c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
         CSTTT(I,J)=1.
      ENDDO
      ENDDO
C     ******************************************************************
c
      TS11=TGBAT0
      QS11=3.799E3/PSFC*EXP(17.26939-4098.026247/(TS11-35.86))*.925
      write(6,66666) ts11,qs11
      IF (ISFC .EQ. 1 .AND. ITOGA .EQ. 1) THEN
c$doacross local(j,i)
          do j=1,jmax
          do i=1,imax
             XYP(i,j)=uu1(i,j,2)
             XY(i,j)=vv1(i,j,2)
             XYM(i,j)=dpt(i,j,2)
             XYPM(i,j)=dqv(i,j,2)
          enddo
          enddo
         call pblin (xyp,xy,xym,xypm,ri180)
      ENDIF
c
C     ******************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         tsdt=0.
         qsdt=0.
         thsdt=0.
         psdt=0.
         do k=1,nz
           tls(k)=0.0
           qls(k)=0.0
           UBT(K)=0.
           VBT(K)=0.
           WBT(K)=0.
           Q1T(K)=0.
           Q2T(K)=0.
           Q1ZT(K)=0.
           Q2ZT(K)=0.
           aq1t(k)=q1_6h(k,1)
           aq2t(k)=q2_6h(k,1)
           aq1zt(k)=q1z_6h(k,1)
           aq2zt(k)=q2z_6h(k,1)
         enddo
         IF (ilif .EQ. 0) Then
           do k=1,nz
             UBT(K)=(UB_6H(K,2)-UB_6H(K,1))*R36
             VBT(K)=(VB_6H(K,2)-VB_6H(K,1))*R36
             WBT(K)=(WB_6H(K,2)-WB_6H(K,1))*R36
             Q1T(K)=(Q1_6H(K,2)-Q1_6H(K,1))*R36
             Q2T(K)=(Q2_6H(K,2)-Q2_6H(K,1))*R36
             Q1ZT(K)=(Q1Z_6H(K,2)-Q1Z_6H(K,1))*R36
             Q2ZT(K)=(Q2Z_6H(K,2)-Q2Z_6H(K,1))*R36
c            aq1t(k)=0.0
c            aq2t(k)=0.0
c            aq1zt(k)=0.0
c            aq2zt(k)=0.0
             if (incar .eq. 1) then
               Q1ZT(K)=0.0
               Q2ZT(K)=0.0
             endif
           enddo
         ENDIF
         tsdt=(ts(2)-ts(1))*R36
         qsdt=(qs(2)-qs(1))*R36
         thsdt=(ths(2)-ths(1))*R36
         psdt=(pss(2)-pss(1))*R36
c$doacross local(j,i)
         do j=1,jmax
           do i=1,imax
             tairsfc(i,j)=ts(1)
             qairsfc(i,j)=qs(1)
             thairsf(i,j)=ths(1)
             pairsfc(i,j)=pss(1)
           enddo
         enddo
cccccccccccccccc
c      istrt_2h=0
      ir2hour=3600
      nudge=0
      r2hour=1./float(ir2hour)
      factor_nuding=r2hour
C     ********
       do k=1,nz
         ubi(k)=0.0
         vbi(k)=0.0
       enddo
c$doacross local(k,i)
       do k=1,itt3
          do i=1,nz
             ub_2h(i,k)=0.
             vb_2h(i,k)=0.
          enddo
       enddo

       do ii=1,itt
          ip=ii+1
          if(ii .eq. itt) ip=itt
          iip3=(ii-1)*3+1
          iipp3=iip3+1
          iippp3=iipp3+1
         do k=1,nz
           y1(k)=(ub_6h(k,ip)-ub_6h(k,ii))/3.
           y2(k)=(vb_6h(k,ip)-vb_6h(k,ii))/3.
           ub_2h(k,iip3)=ub_6h(k,ii)
           ub_2h(k,iipp3)=ub_6h(k,ii)+y1(k)
           ub_2h(k,iippp3)=ub_6h(k,ii)+2.*y1(k)
           vb_2h(k,iip3)=vb_6h(k,ii)
           vb_2h(k,iipp3)=vb_6h(k,ii)+y2(k)
           vb_2h(k,iippp3)=vb_6h(k,ii)+2.*y2(k)
         enddo
       enddo

ctao (12-23-98)
        do k=1,nz
          ubi(k)=ub(k)
          vbi(k)=vb(k)
c         ubi(k)=ub_2h(k,2)
c         vbi(k)=vb_2h(k,2)
        enddo
ctao (12-23-98)
       do ii=1,itt
          iip3=(ii-1)*3+1
          iipp3=iip3+1
          iippp3=iipp3+1
         print*
         print*,'itt =',ii,iip3,iipp3,iippp3
         print*
         write(6,24100)
        do k1=2,kles
          k=kles+2-k1
          km=k-1
          write(6,20300) km,ub_2h(k,iip3),ub_2h(k,iipp3),ub_2h(k,iippp3)
     1                  ,vb_2h(k,iip3),vb_2h(k,iipp3),vb_2h(k,iippp3)
         enddo
       enddo

20300 format(i4,10f12.5,/)
24100 format(//,'level          u           u1          u2          v 
     1       v1          v2')


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         rm0=0.5
         itrn=8*60*60
         rm1=0.25
         itrn1=4*24*60*60
         rm2=0.25
         itrn2=4*24*60*60+12*60*60
         rm3=0.25
          nrun = 2
         ised1=ised(nrun)
       NRAN=0
       RM=RM0/150.
      DO J=JRMS,JRMF
      DO I=IRMS,IRMF
         RX(I,J)=(RAN3(ised1)-.5)*rm
      ENDDO
      ENDDO

        icheck=0

       NRAN=NRAN+1
C     ******************************************************************
      IF (IRS.EQ.0 .AND. N.EQ.1) GO TO 6666
      IF (IRS.EQ.0 .AND. N.EQ.0) GO TO 5555
C
       CALL RESTART ('READ')

C
       WRITE(6,12345) N,ISEC,NRAN,SEC,AMINUT,RDT
        DT=RDT
        D2T=DT+DT
C
       CALL SEPCA
C
       rm=rm0/150.
        if(isec.ge.itrn) rm=rm1/150.
        if(isec.ge.itrn1) rm=rm2/150.
        if(isec.ge.itrn2) rm=rm3/150.

      DO IJKR=1,NRAN-1
         DO J=JRMS,JRMF
            DO I=IRMS,IRMF
               RX(I,J)=(RAN3(sed1)-.5)*rm
            ENDDO
         ENDDO
      ENDDO

C     ******************************************************************
 3333 N=N+1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      rbud = rbud + 1.

      call timing (time2)
      dtime=time2-time1
      time1=time2

      if(isec.eq.isec/dt*dt) then
        AMINT=AMINUT-DT/60.
        AHRT=AMINT/60.
        PRINT 241, N-1,SEC-DT,AMINT,AHRT,DTIME
  241   FORMAT('TIMESTEP',I5,' MODEL TIME(S,M,H)=',F10.2,
     $       F10.2,F7.2,
     $       3X,'CPU(S)=',F8.3)
      endif

       cpu=cpu+dtime

      SEC=SEC+DT
      ISEC=SEC+.1
      ISEC1=ISEC
      AMINUT=SEC/60.
cc      ISEC1=ISEC
      ID=0
       IF (ISEC .EQ. ISEC/ISTATMIN*ISTATMIN) ID=1
       IF (ITEST .EQ. 0) ID=0
c
      idpbl=0
        if(isec.eq.isec/IPBLDT*IPBLDT) idpbl=1
        if (isfc .eq. 0) idpbl=0
        if (itoga .ne. 1) idpbl=0
c
       HRL=HRL+dt/3600.
      IRADD=0
      IF (IRAD .EQ. 1) THEN
        IF(ISEC.EQ.ISEC/IRADDT*IRADDT) THEN
          call zang (cosz,rlat,rmonth,riday,hrl)
          iradd=1
        ENDIF
      ENDIF
c
        if (icase .eq. 0) then
          ICY_UB=ISEC/ir36+1
          DO ICY=1,ICY_UB
            ICHK=ISTRT+(ICY-1)*ir36
            IF(ISEC.EQ.ICHK)THEN
              ITFOR=ICY
c             IF (ilif .EQ. 0) ITFOR=ICY+1
              IF (ilif .EQ. 0)ITFOR=MIN(ICY+1,ITT-1) ! shie 11/16/01 3d
              IF (IWBAR .EQ. 0) THEN
                DO K=2,KLES
                  UBT(K)=(UB_6H(K,itfor+1)-UB_6H(K,itfor))*R36
                  VBT(K)=(VB_6H(K,itfor+1)-VB_6H(K,itfor))*R36
                  WBT(K)=(WB_6H(K,itfor+1)-WB_6H(K,itfor))*R36
                  Q1T(K)=(Q1_6H(K,itfor+1)-Q1_6H(K,itfor))*R36
                  Q2T(K)=(Q2_6H(K,itfor+1)-Q2_6H(K,itfor))*R36
                  Q1ZT(K)=(Q1Z_6H(K,itfor+1)-Q1Z_6H(K,itfor))*R36
                  Q2ZT(K)=(Q2Z_6H(K,itfor+1)-Q2Z_6H(K,itfor))*R36
                  aq1t(k)=q1_6h(k,itfor)
                  aq2t(k)=q2_6h(k,itfor)
                  aq1zt(k)=q1z_6h(k,itfor)
                  aq2zt(k)=q2z_6h(k,itfor)
                ENDDO
                tsdt=(ts(itfor+1)-ts(itfor))*R36
                qsdt=(qs(itfor+1)-qs(itfor))*R36
                thsdt=(ths(itfor+1)-ths(itfor))*R36
                psdt=(pss(itfor+1)-pss(itfor))*R36
c$doacross local(j,i)
                do j=1,jmax
                  do i=1,imax
                    tairsfc(i,j)=ts(itfor)
                    qairsfc(i,j)=qs(itfor)
                    thairsf(i,j)=ths(itfor)
                    pairsfc(i,j)=pss(itfor)
                  enddo
                enddo
              ELSE
                DO K=2,KLES
                  if (iuvbar .eq. 0) UBT(K)=0.
                  if (iuvbar .eq. 0) VBT(K)=0.
                  WBT(K)=0.
                  Q1T(K)=0.
                  Q2T(K)=0.
                ENDDO
                tsdt=0.
                qsdt=0.
                thsdt=0.
                psdt=0.
              ENDIF
              WRITE (6,*) 'UPDATE UBT AT ', ISEC,' ICY=',itfor
              WRITE (6,*) 'UPDATE WBT AT ', ISEC,' ICY=',itfor
            ENDIF
          ENDDO
        endif
C-----------------------------------------------------------------------
ccc
      IF (IDPBL .EQ. 1) THEN

          if (iricont .eq. 0) then
c$doacross local(j,i)
          do j=2,jles
          do i=2,iles
            riold180(i,j) = 0.
          enddo
          enddo

          else

c$doacross local(j,i)
          do j=2,jles
          do i=2,iles
            riold180(i,j)=ri180(i,j)/float(iricont)
          enddo
          enddo

          endif

c$doacross local(j,i)
          do j=1,jmax
          do i=1,imax
             XYP(i,j)=uu1(i,j,2)
             XY(i,j)=vv1(i,j,2)
             XYM(i,j)=dpt(i,j,2)
             XYPM(i,j)=dqv(i,j,2)
          enddo
          enddo
c
          IF(IJKADV .EQ. 1) THEN
c$doacross local(j,i)
             do j=1,jmax
             do i=1,imax
                xyp(i,j)=.5*(3.*uu1(i,j,2)-umd(i,j,2))
                xy(i,j)=.5*(3.*vv1(i,j,2)-vmd(i,j,2))
             enddo
             enddo
          endif

         call pblin (xyp,xy,xym,xypm,ri180)

          iricont=0
c$doacross local(j,i)
          do j=2,jles
          do i=2,iles
             ri180(i,j)=0.
          enddo
          enddo

      ENDIF
c
cc
ccc
C*****************************************************************************
      TA11=TA1(1)
      QA11=QA1(1)
      TA12=TA1(2)
      QA12=QA1(2)
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
       XYM(I,J)=(DPT1(I,J,1)+TA11)*(1.+.61*(DQV1(I,J,1)+QA11)
     1   -QCL1(I,J,1)-QRN1(I,J,1)-QCI1(I,J,1)-QCS1(I,J,1)-QCG1(I,J,1))

       XY(I,J)=(DPT1(I,J,2)+TA12)*(1.+.61*(DQV1(I,J,2)+QA12)
     1   -QCL1(I,J,2)-QRN1(I,J,2)-QCI1(I,J,2)-QCS1(I,J,2)-QCG1(I,J,2))
      ENDDO
      ENDDO

      DO K=2,KLES
       KP=K+1
       TA1KP=TA1(KP)
       QA1KP=QA1(KP)
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
       XYP(I,J)=(DPT1(I,J,KP)+TA1KP)*(1.+.61*(DQV1(I,J,KP)+QA1KP)
     1  -QCL1(I,J,KP)-QRN1(I,J,KP)-QCI1(I,J,KP)-QCS1(I,J,KP)
     2  -QCG1(I,J,KP))
       IF(XYP(I,J)-XYM(I,J) .GE. 0.0) THEN
        XXW(I,J,K)=0.0
       ELSE
        XXW(I,J,K)=0.5
       ENDIF
       IF(K.NE.KLES) THEN
        XYM(I,J)=XY(I,J)
        XY(I,J)=XYP(I,J)
       ENDIF
      ENDDO
      ENDDO
      ENDDO
c new
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
       DF0(I,J,2)=0.0
       DFC(I,J,2)=0.0
       XXW(I,J,1)=XXW(I,J,2)
       XXW(I,J,KMAX)=XXW(I,J,KLES)
      ENDDO
      ENDDO

      DO K=2,KMAX
       Y1(K)=F0(K)+F0(K-1)
       Y2(K)=.5*Y1(K)
      ENDDO

      DO K=3,KMAX
       KM=K-1
       Y1K=2.17375*Y1(K)
       QAK=QA(K)+QA(K-1)
       TAK=TA(K)+TA(K-1)
       TA1K=TA1(K)-TA1(K-1)
       QA1K=QA1(K)-QA1(K-1)
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
        DF0(I,J,K)=Y1K*(DQV(I,J,K)+DQV(I,J,KM)+QAK)
     2      /(DPT(I,J,K)+DPT(I,J,KM)+TAK)**2

         DFC(I,J,K)=-AM1(K)*(AK(I,J,K)*(1.+XXW(I,J,K))+AK(I,J,KM)
     1     *(1.+XXW(I,J,KM)))*(DPT1(I,J,K)-DPT1(I,J,KM)+TA1K
     2     +Y2(K)*(DQV1(I,J,K)-DQV1(I,J,KM)+QA1K))
     3     /((1.+DF0(I,J,K)*Y2(K))*DZ)
       ENDDO
       ENDDO
      ENDDO

c      DO K=3,KMAX
c        KM=K-1
c         QAK=QA(K)+QA(K-1)
c         TAK=TA(K)+TA(K-1)
c         TA1K=TA1(K)-TA1(K-1)
c         QA1K=QA1(K)-QA1(K-1)
c         Y1K=2.17375*Y1(K)
c       DO J=2,JLES
c       DO I=2,ILES
c         DF0(I,J,K)=Y1K
c     1      *(DQV(I,J,K)+DQV(I,J,KM)+QAK)
c     2      /(DPT(I,J,K)+DPT(I,J,KM)+TAK)**2
c         DFC(I,J,K)=-AM1(K)*(AK(I,J,K)*(1.+XXW(I,J,K))+AK(I,J,KM)
c     1     *(1.+XXW(I,J,KM)))*(DPT1(I,J,K)-DPT1(I,J,KM)+TA1K
c     2     +Y2(K)*(DQV1(I,J,K)-DQV1(I,J,KM)+QA1K))
c     3     /((1.+DF0(I,J,K)*Y2(K))*DZ)
c       ENDDO
c       ENDDO
c      ENDDO


C     *******************************************************
      ISMG=1
      kt1=2
      kt2=kles
      irf=1
      irfg=0
      idq=0
      IADVH=IADVH1
        IF (IJKADV .EQ. 1)  IADVH=2

      icheck=0
      CALL ADVECT (DPT,DPT1,TA1,ID,ID, 0, 0, ID, 1., 1.)
C     **********************************************************

c$doacross local(k,j,i)
      DO K=3,KMAX
      DO J=2,JLES
      DO I=2,ILES
         DFC(I,J,K)=DFC(I,J,K)*DF0(I,J,K)
      ENDDO
      ENDDO
      ENDDO


C     **********************************************************
      ISMG=2
      icheck=1
      CALL ADVECT (DQV,DQV1,QA1, 0, 0, ID, 0, 0, 1., 1.)
C     **********************************************************

      DO K=3,KMAX
        KM=K-1
        A0K=AM1(K)*RDZ
        QA1K=QA1(K)
        QA1KM=QA1(KM)
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
       DFC(I,J,K)=-DFC(I,J,K)-A0K*(AK(I,J,K)+AK(I,J,KM))
     1            *(DQV1(I,J,K)-DQV1(I,J,KM)+QA1K-QA1KM
     2      +QCL1(I,J,K)-QCL1(I,J,KM)+QCI1(I,J,K)-QCI1(I,J,KM))
      ENDDO
      ENDDO
      ENDDO

      IRF=0
      ISMG=3
      idq=id

c$doacross local(k,j,i)
       DO K=1,KMAX
       DO J=1,JMAX
       DO I=1,IMAX
          QCL2(I,J,K)=QCL(I,J,K)
       ENDDO
       ENDDO
       ENDDO

C     *************************************************************
      CALL ADVECT (QCL,QCL1,QAA1, 0, 0, 0, 0, 0, 1., .0)
C     *************************************************************
      IF (ICE .EQ. 1) THEN
         ISMG=4
         CALL ADVECT (QCI,QCI1,QAA1, 0, 0, 0, 0, 0, .0, 1.)
      ENDIF
C     *************************************************************



c$doacross local(k,j,i)
       DO K=1,KMAX
       DO J=1,JMAX
       DO I=1,IMAX
          XXX(I,J,K)=WW1(I,J,K)
       ENDDO
       ENDDO
       ENDDO

      IRSG=0
      IF (ICE .EQ. 1) THEN
         IF (ICERH.EQ.1) CALL TERVRH (IRSG,RHO1,FV1)
      ELSE
      ENDIF
       A1=36.E3*RHO1(2)
       A2=2.77778E-4*DT

c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
         RI(I,J)=A1*QRN(I,J,2)*DF0(I,J,2)
         ri180(i,j)=ri180(i,j)+ri(i,j)
         AR(I,J)=AR(I,J)+A2*RI(I,J)
      ENDDO
      ENDDO

      IRICONT=IRICONT+1
      IRICONT1=IRICONT1+1
C
      IF (IJKADV .EQ. 0) THEN
c$doacross local(k,j,i)
         DO K=1,KMAX
         DO J=1,JMAX
         DO I=1,IMAX
            WW1(I,J,K)=XXX(I,J,K)-DF0(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      ISMG=5
      CALL ADVECT (QRN,QRN1,QAA1, 0, 0, 0, 1, 0, 0., 0.)
C     ************************************************************
      IF (ICE .EQ. 1) THEN
         IRSG=1
         IF (ICERH.EQ.1) CALL TERVRH (IRSG,RHO1,FV1)
         IF (IJKADV .EQ. 0) THEN
c$doacross local (k,j,i)
            DO K=1,KMAX
            DO J=1,JMAX
            DO I=1,IMAX
               WW1(I,J,K)=XXX(I,J,K)-DF0(I,J,K)
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         CALL ADVECT (QCS,QCS1,QAA1, 0, 0, 0, 1, 0, 0., 0.)
C     *************************************************************
         IRSG=2
         IF (ICERH.EQ.1) CALL TERVRH (IRSG,RHO1,FV1)
         IF (IJKADV .EQ. 0) THEN
c$doacross local(k,j,i)
            DO K=2,KMAX
            DO J=2,JLES
            DO I=2,ILES
               WW1(I,J,K)=XXX(I,J,K)-DF0(I,J,K)
            ENDDO
            ENDDO
            ENDDO
          ENDIF
         CALL ADVECT (QCG,QCG1,QAA1, 0, 0, 0, 1, 0, 0., 0.)
      ENDIF
       icheck=0
C     *************************************************************


c$doacross local(k,j,i)
       DO K=1,KMAX
       DO J=1,JMAX
       DO I=1,IMAX
          WW1(I,J,K)=XXX(I,J,K)
       ENDDO
       ENDDO
       ENDDO

      IRF=1
      SMR=0.

c$doacross local(j,i),reduction(smr)
      DO J=2,JLES
      DO I=2,ILES
         SMR=SMR+AR(I,J)
      ENDDO
      ENDDO

      SMR=SMR*RIJL2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       D22T=D2T
         IF(IJKADV .EQ. 1) D22T=DT
        do k=2,kles
          y1(k)=d22t*rdz*rrho(k)
          y2(k)=.5*d22t*fcorg*rdz*am(k)*tb(k)
          IF (ICASE .EQ. 1) THEN
              Q1ZT(K)=tlsw(k)*1.004*3600.0*pi(k)
              Q2ZT(K)=qlsw(k)*2500.0*3600.0
              Q1T(K)=0.
              Q2T(K)=0.
              if (isui .eq. 1) Q1ZT(K)=0.
              if (isui .eq. 1) Q2ZT(K)=0.
          ENDIF
        enddo
        rtfor=d22t/(24.*3600.)
        rqfor=d22t/(24.*3600.)
         if (incar .eq. 0) then
            rqfor=d22t/(2500.*24.*3600.)
         endif
        afact=1.
        if (ijkadv .eq. 0) afact=0.5


        do k=2,kles
           kp=k+1
           km=k-1
           a3=0.
           a4=0.
           a5=0.
           a6=0.
           a7=0.
           a8=0.
           am1fd=am1(k)*fd(k)
           am1fdp=am1(kp)*fd(kp)
           am1fe=am1(k)*fe(k)
           am1fep=am1(kp)*fe(kp)
           y1k=y1(k)
           y2k=y2(k)
           pik=1./pi(k)
           aq1tk=aq1t(k)*rtfor*pik
           aq2tk=aq2t(k)*rqfor
           aq1ztk=aq1zt(k)*rtfor*pik
           aq2ztk=aq2zt(k)*rqfor
        IF (ISNGTAO .EQ. 0) THEN

           DO J=2,JLES
           DO I=2,ILES
              IF(IJKADV .EQ. 1) THEN
                 y3i=.5*(3.*ww1(i,j,k)-wmd(i,j,k))
                 y4i=.5*(3.*ww1(i,j,kp)-wmd(i,j,kp))
              ELSE
                 Y3I=WW1(I,j,K)
                 Y4I=WW1(I,j,KP)
              ENDIF

          DPT(I,j,K)=DPT(I,j,K)-(AM1FD*Y3I+AM1FDP*Y4I)*Y1K-y0(i,j)*aQ1TK
          DQV(I,j,K)=DQV(I,j,K)-(AM1FE*Y3I+AM1FEP*Y4I)*Y1K-y0(i,j)*aQ2TK

          a77=-y0(i,j)*aq1ztk
          a88=-y0(i,j)*aq2ztk
          a1=-(am1fd*y3i+am1fdp*y4i)*y1k-y0(i,j)*aq1tk-a77
          a2=-(am1fe*y3i+am1fep*y4i)*y1k-y0(i,j)*aq2tk-a88


c             q1_g_v(i,j,k)=q1_g_v(i,j,k)-(AM1FD*Y3I
c     1                    +AM1FDP*Y4I)*Y1K
c
c             q2_g_v(i,j,k)=q2_g_v(i,j,k)-(AM1FE*Y3I
c     1                    +AM1FEP*Y4I)*Y1K
c              DPT(I,j,K)=DPT(I,j,K)-(AM1FD*Y3I+AM1FDP*Y4I)*Y1K
c    2              +Y2K*(VV1(I,K)-VB1(K))*(UU1(I,j,KP)-UU1(I,j,KM))
c    3              -y0(i,j)*aQ1TK-y0(i,j)*aQ1ZTK
c
c      DQV(I,j,K)=DQV(I,j,K)-(AM1FE*Y3I+AM1FEP*Y4I)*Y1K-y0(i,j)*aQ2TK
c    3              -y0(i,j)*aQ2ZTK

          a5=a5-(am1fd*y3i+am1fdp*y4i)*y1k
          a6=a6-(am1fe*y3i+am1fep*y4i)*y1k
          a7=a7+a77
          a8=a8+a88
          a3=a3+a1
          a4=a4+a2

          ENDDO
          ENDDO


          ELSE

         DO J=2,JLES
         DO I=2,ILES
          IF(IJKADV .EQ. 1) THEN
            y3(i)=.5*(3.*ww1(i,j,k)-wmd(i,j,k))
            y4(i)=.5*(3.*ww1(i,j,kp)-wmd(i,j,kp))
          ELSE
            Y3(I)=WW1(I,j,K)
            Y4(I)=WW1(I,j,KP)
          ENDIF
            dpt(i,j,k)=dpt(i,j,k)-(am1fd*(Y3(I)+wb(k)*y0(i,j))
     1                +am1fdp*(Y4(I)+wb(kp)*y0(i,j)))*y1k
c    2                +y2(k)*(vv1(i,k)-vb1(k))*(uu1(i,kp)-uu1(i,km))
     3                -y0(i,j)*aQ1TK
            dqv(i,j,k)=dqv(i,j,k)-(am1fe*(Y3(I)+wb(k)*y0(i,j))
     1                +am1fep*(Y4(I)+wb(kp)*y0(i,j)))*y1k
     2                -y0(i,j)*aQ2TK
             a1=-(am1fd*(Y3(I)+wb(k)*y0(i,j))
     1          +am1fdp*(Y4(I)+wb(kp)*y0(i,j)))*y1k
     2          -y0(i,j)*aQ1TK
             a2=-(am1fe*(Y3(I)+wb(k)*y0(i,j))
     1          +am1fep*(Y4(I)+wb(kp)*y0(i,j)))*y1k
     2          -y0(i,j)*aQ2TK
             a77=-(am1fd*wb(k)*y0(i,j)+am1fdp*wb(kp)*y0(i,j))*y1k
             a88=-(am1fe*wb(k)*y0(i,j)+am1fep*wb(kp)*y0(i,j))*y1k

          a5=a5-(am1fd*y3(i)+am1fdp*y4(i))*y1k
          a6=a6-(am1fe*y3(i)+am1fep*y4(i))*y1k
          a7=a7+a77
          a8=a8+a88
          a3=a3+a1
          a4=a4+a2
        ENDDO
        ENDDO

          ENDIF

        tlsw(k)=tlsw(k)+a5*rijl2*afact
        qlsw(k)=qlsw(k)+a6*rijl2*afact
        ttlsw(k)=ttlsw(k)+a7*rijl2*afact
        tqlsw(k)=tqlsw(k)+a8*rijl2*afact
        tls(k)=tls(k)+a3*rijl2*afact
        qls(k)=qls(k)+a4*rijl2*afact
      ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (ISEC.ne.ISEC/300*300) GO TO 6666
      IF (ISEC.GE.ITRN) RM=RM1/150.
      if (isec.ge.itrn1) rm=rm2/150.
      if (isec.ge.itrn2) rm=rm3/150.
       DO J=JRMS,JRMF
       DO I=IRMS,IRMF
          XY(I,J)=RX(I,J)
c         RX(I,J)=(RANF()-.5)*RM
          RX(I,J)=(RAN3(ised1)-.5)*RM
       ENDDO
       ENDDO
       IF (IJKADV .EQ. 0) THEN
          DO J=JRMS,JRMF
          DO I=IRMS,IRMF
            DPT(I,J,kdist)=DPT(I,J,kdist)+(XY(I,J)-RX(I,J))*DT
          ENDDO
          ENDDO
       ENDIF
       NRAN=NRAN+1
C     *************************************************************
 6666 CONTINUE
       D22T=D2T
         IF(IJKADV .EQ. 1) D22T=DT

c$doacross local(j,i)
      DO J=JRMS,JRMF
      DO I=IRMS,IRMF
         DPT(I,J,kdist)=DPT(I,J,kdist)+RX(I,J)*D22T
      ENDDO
      ENDDO
C
CC    ***   PARAMERTERIZED SURFACE FLUXES   ****************************

      IF (itoga .eq. 2) THEN
        D22T=D2T
          IF(IJKADV .EQ. 1) D22T=DT
       A1=1.E-3*RHO1(2)
       A2=D22T*AM(2)*RDZ*RRHO(2)
       A3=.0
       A4=.0

        DO J=2,JLES
          JP=J+1
          IF (J .EQ. JLES) JP=2
        DO I=2,iles
            IP=I+1
            IF (I .EQ. ILES) IP=2
cccccc
        IF (isfcave .EQ. 0) THEN
          VVSS=SQRT( (.5*(UU1(I,J,2)+UU1(I+1,J,2))+UBAR)**2 
     1             + (.5*(VV1(I,J,2)+VV1(I,Jp,2))+VBAR)**2 )
          IF(IJKADV .EQ. 1) VVSS=SQRT((.25*(3.*UU1(I,J,2)-UMD(I,J,2)
     1                            +3.*UU1(IP,J,2)-UMD(IP,J,2))+ubar)**2
     2                            + (.25*(3.*VV1(I,J,2)-VMD(I,J,2)
     3                            +3.*VV1(I,JP,2)-VMD(I,JP,2))+vbar)**2)
        ELSE
          VVSS=SQRT( (UB1(2)+ubar)**2 + (VB1(2)+vbar)**2 )
        ENDIF
cccccc
            Y4(I)=MAX(VVSS, 100.)
           Y1(I)=A1*(1.1+.04*(.01*Y4(I)))*Y4(I)
           if (itest_model .eq. 1) then
             Y2(I)=MAX(Y1(I)*(TS11-TA1(2)-DPT(I,J,2)), 0.)
             Y3(I)=MAX(Y1(I)*(QS11-QA1(2)-DQV(I,J,2)), 0.)
           else
             Y2(I)=MAX(Y1(I)*(thairsf(I,j)-TA1(2)-DPT(I,J,2)), 0.)
             Y3(I)=MAX(Y1(I)*(qairsfc(I,j)-QA1(2)-DQV(I,J,2)), 0.)
             IF (isfcave .EQ. 0) THEN
                y2(i)=max(y1(i)*(thairsf(i,j)-ta1(2)) ,0.)
                y3(i)=max(y1(i)*(qairsfc(i,j)-qa1(2)) ,0.)
             ENDIF
           endif

c         q1_d_v(i,j,2)=q1_d_v(i,j,2)+a2*y2(i)
c         q2_d_v(i,j,2)=q2_d_v(i,j,2)+a2*y3(i)

         DPT(I,J,2)=DPT(I,J,2)+A2*Y2(I)
         DQV(I,J,2)=DQV(I,J,2)+A2*Y3(I)
          A3=A3+Y2(I)
          A4=A4+Y3(I)
C SAVE SENSIBLE AND LATENT HEAT FLUXES
         TS0(I,J)=A2*Y2(I)/DT
         QSS0(I,J)=A2*Y3(I)/DT
        SUN_4(I,j,3)=Y2(I)
        SUN_4(I,j,4)=Y3(I)
         tsfc_1(i,j)=y2(i)
         qsfc_1(i,j)=y3(i)
        ENDDO
        ENDDO
        STF(2)=STF(2)+A3*DT*rijl2
        SQF(2)=SQF(2)+A4*DT*rijl2
      ENDIF

CC    ***   INSERT RADIATIVE TRANSFER PROCESSES HERE *****
      IF (IRAD.EQ.1) THEN
c         iflag=1
         ircl=ircl+1
c         if(isec.eq.isec/10800*10800) iflag=2
         IF(IRADD.EQ.1)  CALL RADRAT (1)
c         if (iradd .eq. 1) call radrat (iflag,cosz,npp1,pl,pa)
         D22T=D2T
            IF(IJKADV .EQ. 1) D22T=DT
c$doacross local(k,j,i)
         DO K=2,KLES
         DO J=2,JLES
         DO I=2,ILES
c            q1_rad(i,j,k)=q1_rad(i,j,k)+(RSW(I,J,K)-RLW(I,J,K))*D22T
            DPT(I,J,K)=DPT(I,J,K)+(RSW(I,J,K)-RLW(I,J,K))*D22T
         ENDDO
         ENDDO
         ENDDO
      ENDIF
CCCCCC   *************************************************
       ISEC1=ISEC+DT
      IF(IWATER .EQ. 1) THEN
        CALL DOMAINW (1)
        CALL SATICEW ( )
      ELSE
        CALL DOMAIN (1)
        IF (ICERH .EQ. 1) CALL SATICERH (FV)

        IF (ICELIN .EQ. 1) CALL SATICEL (FV)

      ENDIF
CC   ***  CHURCHILL AND HOUZES' CONVECTIVE AND ANVIL SEPARATION   *****
      IF (ISEC .EQ. ISEC/ISTATMIN*ISTATMIN) CALL SEPCA 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do k=2,kles
          wb(k)=wb(k)+wbt(k)*dt
ctao (12-23-98)
          ubi(k)=ubi(k)+ubt(k)*dt
          vbi(k)=vbi(k)+vbt(k)*dt
ctao (12-23-98)
          aq1t(k)=aq1t(k)+q1t(k)*dt
          aq2t(k)=aq2t(k)+q2t(k)*dt
          aq1zt(k)=aq1zt(k)+q1zt(k)*dt
          aq2zt(k)=aq2zt(k)+q2zt(k)*dt
          if (incar .eq. 1) then
             aq1zt(k)=0.0
             aq2zt(k)=0.0
          endif
        enddo
c$doacross local(j,i)
        do j=1,jmax
        do i=1,imax
          tairsfc(i,j)=tairsfc(i,j)+tsdt*dt
          qairsfc(i,j)=qairsfc(i,j)+qsdt*dt
          thairsf(i,j)=thairsf(i,j)+thsdt*dt
          pairsfc(i,j)=pairsfc(i,j)+psdt*dt
        enddo
        enddo
          iihalf=il2/2+1
          jjhalf=jl2/2+1
         if(isec.eq. isec/10800*10800)
     1       print*,aminut,tairsfc(iihalf,jjhalf),
     2              qairsfc(iihalf,jjhalf)*1000.,
     2              pairsfc(iihalf,jjhalf)/1000.,thairsf(iihalf,jjhalf)

      D2T=DT+DT
      ISEC=INT(SEC+.1)
      NSML=INT(D2T/DTS)
      ISEC1=ISEC
      SMT=0.0
      SMV=0.0
      SMQ=0.0
      SMC=0.0
      SME=0.0
      SMQC=0.0
      SMQR=0.0
      SMQI=0.0
      SMQS=0.0
      SMQG=0.0
c
      sddd=0.
      ssss=0.
      shhh=0.
      sccc=0.
      smmm=0.
      sfff=0.
c
      DO K=2,KLES
        Y1(K)=10.*DZ*RHO(K)/AM(K)
        Y2(K)=Y1(K)/F0(K)
c
       sddd=sddd+s_dep(k)*y1(k)
       ssss=ssss+s_sub(k)*y1(k)
       shhh=shhh+s_qrs(k)*y1(k)
       sccc=sccc+s_qrl(k)*y1(k)
       smmm=smmm+s_mel(k)*y1(k)
       sfff=sfff+s_frz(k)*y1(k)
c
       SMT=SMT+ST(K)*Y2(K)
       SMV=SMV+SV(K)*Y1(K)
       SMQ=SMQ+SQ(K)*Y1(K)
       SMC=SMC+SC(K)*Y1(K)
       SME=SME+SE(K)*Y1(K)
       SMQC=SMQC+SQC(K)*Y1(K)
       SMQR=SMQR+SQR(K)*Y1(K)
       SMQI=SMQI+SQI(K)*Y1(K)
       SMQS=SMQS+SQS(K)*Y1(K)
       SMQG=SMQG+SQG(K)*Y1(K)
       ENDDO
C     ******      *********************
      IF (ISEC.GE.IWTP) THEN
         IF(ISEC.EQ.ISEC/IMICRO*IMICRO) CALL WMICRO(ifmicro)
         IF(ISEC.EQ.ISEC/IDYN*IDYN) CALL WDYN(ifdyn)
         IF(ISEC.EQ.ISEC/ISTAT*ISTAT) CALL WSTAT(ifstat)
      ENDIF
C     ******      *********************
      IF (ISEC.EQ.ISEC/IMP*IMP) THEN
         icot=1
         jcot=1
         if (nxi .eq. 128) icot=1
         if (nxi .eq. 256) icot=2
         if (nxi .eq. 512) icot=4
         if (nxi .eq. 1024) icot=8
         if (nyj .eq. 128) jcot=1
         if (nyj .eq. 256) jcot=2
         if (nyj .eq. 512) jcot=4
         if (nyj .eq. 1024) jcot=8
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL WCLOUD (N,AMINUT,icot,jcot,CCON,SMT,SMV,SMQ,
     1                SMC,SME,SMR,SMQC,SMQR,SMQI,SMQS,SMQG,
     2                sddd,ssss,shhh,sccc,smmm,sfff)
         CALL WCLOUDW (N,AMINUT,icot,jcot)
      ENDIF
c      *****************************************************************
       CALL WTAPRI(ifrain)
       IF (ISEC.EQ.ISEC/IMX*IMX) THEN
        CALL TMAX (WW1,WP,AMINUT,WMAX,WMIN,1)
        CALL TMAX (uu1,uP,AMINUT,WMAX,WMIN,1)
        CALL TMAX (vv1,vP,AMINUT,WMAX,WMIN,1)
         M=M+1
        RWMAX(M)=WMAX
        RWMIN(M)=WMIN
        CALL TMAX (DPT,PT,AMINUT,WMAX,WMIN,2)
        CALL TMAX (DQV,SH,AMINUT,WMAX,WMIN,3)

c$doacross local(k,qbk,j,i)
        DO K=1,KMAX
        QBK=QB(K)
        DO J=1,JMAX
        DO I=1,IMAX
           DF0(I,J,K)=QBK+DQV(I,J,K)
        ENDDO
        ENDDO
        ENDDO

c       CALL TMAX (DF0,SH,AMINUT,WMAX,WMIN,3)
        CALL TMAX (DF0,SHd,AMINUT,WMAX,WMIN,3)  ! cccshie 11/16/01
        CALL TMAX (QCL,CD,AMINUT,WMAX,WMIN,3)
        CALL TMAX (QRN,RN,AMINUT,WMAX,WMIN,3)
        if (iwater .eq. 0) then
          CALL TMAX (QCI,QI,AMINUT,WMAX,WMIN,3)
          CALL TMAX (QCS,RS,AMINUT,WMAX,WMIN,3)
          CALL TMAX (QCG,RG,AMINUT,WMAX,WMIN,3)
        endif
        CALL TMAX ( AK,DF,AMINUT,WMAX,WMIN,4)
C       if (isec .eq. isec/600*600) then
C        call page (dpt,pt,aminut)
C        call page (dqv,sh,aminut)
C        call page (ww1,wp,aminut)
C        call page (uu1,up,aminut)
C        call page (vv1,vp,aminut)
C       endif
      ENDIF
C     ***********************************


      IF (ISEC.EQ.ISEC/INP*INP) THEN
      WRITE(6,125)
      DO K1=2,KLES
        K=KLES+2-K1
        WRITE(6,130) K,ST(K),SV(K),SQ(K),SC(K),SE(K),SQA(K),UB(K),VB(K),
     1  WB(K),SQAK(K),SQAAQ(K)
      ENDDO
c
      write(6,5125)
      do k1=2,kles
         k=kles+2-k1
       write(6,134) k,s_dep(k),s_sub(k),s_qrs(k),s_qrl(k),s_mel(k),
     1              s_frz(k)
      enddo
 5125 FORMAT(1X,'LV',6X,'DEP',6X,'SUB',6X,'QRH',6X,'QRC',6X,'MET',6X,'F
     1EZ')
c
      WRITE(6,126)
      DO K1=2,KLES
       K=KLES+2-K1
      WRITE(6,134) K,SQC(K),SQR(K),SQI(K),SQS(K),SQG(K),STQC(K),STQR(K),
     1  STQI(K),STQS(K),STQG(K),wb(k)
      ENDDO
      WRITE(6,138)
      DO K1=2,KLES
       K=KLES+2-K1
      WRITE(6,134) K,SMF(K),SMU(K),SMD(K),SMN(K),CLD(K),CLDU(K),CLDD(K),
     1   CLDN(K),SQF1(K),SQU1(K),SQD1(K),SQN1(K)
      ENDDO
      ENDIF
c     ******************************************************************
      if(isec.eq.isec/21600*21600) then
      write(6,1025)
      do k1=2,kles
        k=kles+2-k1
        km=k-1
       write(6,134) km,Q1ZT(k),Q2ZT(k),Q1T(k),Q2T(k),aq1t(k),aq2t(k)
      ENDDO

       do k=2,kles
         y1(k)=0.
         y2(k)=0.
         y3(k)=0.
         y4(k)=0.
         y5(k)=0.
         y6(k)=0.
       enddo

       do k=2,kles
       do j=2,jles
       do i=2,iles
         y1(k)=y1(k)+ww1(i,j,k)
         y2(k)=y2(k)+uu1(i,j,k)
         y3(k)=y3(k)+vv1(i,j,k)
         y4(k)=y4(k)+w(i,j,k)
         y5(k)=y5(k)+u(i,j,k)
         y6(k)=y6(k)+v(i,j,k)
       enddo
       enddo
       enddo

      do k1=2,kles
       k=kles+2-k1
       km=k-1
        a1=ubt(k)*6.*3600.
        a2=vbt(k)*6.*3600.
        a3=wbt(k)*6.*3600.
        a4=q1t(k)*6.*3600.
        a5=q2t(k)*6.*3600.
        a6=q1zt(k)*6.*3600.
        a7=q2zt(k)*6.*3600.
       write(6,134) km,a1,a2,a3,a4,a5,a6,a7,y1(k)*rijl2,y2(k)*rijl2,
     1              y3(k)*rijl2,y4(k)*rijl2,y5(k)*rijl2,y6(k)*rijl2
      enddo
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 5555 IRF=1
      ISMG=6
       id=0
       idq=0

      IADVH=IADVH1
        IF (IJKADV .EQ. 1) IADVH=2
      CALL AKCOEF

      if (ijkadv .eq. 1) then
c$doacross local(k,j,i)
         DO K=1,KMAX
         DO J=1,JMAX
         DO I=1,IMAX
            UMD(I,J,K)=UU1(I,J,K)
            VMD(I,J,K)=VV1(I,J,K)
            WMD(I,J,K)=WW1(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      endif

      IRF=1
       IADVH=IADVH1
      CALL ADVECTW (SQ)
      CALL ADVECTU (nudge)
      CALL ADVECTV (nudge)
      IRF=0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call cmpf
        CALL SLVPI (1)
C       CALL PRESR (1)
C       CALL POISN (1)

        call cmpwuv

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       IF(ISEC.eq.ISEC/IMX*IMX) THEN
          prefer=dfc(2,2,kles)
        IF (LIPPS .EQ. 1) prefer=dfc(2,2,KLes)/TB(KLes)

          IF (LIPPS .EQ. 0) THEN
c$doacross local(k,j,i)
          do k=1,kmax
          do j=1,jmax
          do i=1,imax
             dfc(i,j,k)=dfc(i,j,k)-prefer
          enddo
          enddo
          enddo
          
          ELSEIF (LIPPS .EQ. 1) THEN
c$doacross local(k,j,i)
          do k=1,kmax
          do j=1,jmax
          do i=1,imax
               dfc(I,j,K)=dfc(I,j,K)/TB(K)-prefer
          enddo
          enddo
          enddo
          ENDIF

         call tmax (dfc,pp,aminut,wmax,wmin,3)
       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC    ***************************
      IF (ISEC.LT.IEND) GO TO 3333
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL WTAPRI(ifrain)
      CALL WTAPS (RWMAX,NT2,ifdyn)
      CALL WTAPS (RWMIN,NT2,ifdyn)
      CALL WTAPS (CNV,NT2,ifdyn)
      CALL WTAPS (ANV,NT2,ifdyn)
      CALL WTAPS (SPN,NT2,ifdyn)
      CALL WSTAT(ifstat)
        RDT=DT
      CALL RESTART ('WRITE')
      WRITE(6,12345) N,ISEC,NRAN,SEC,AMINUT,RDT

      WRITE(6,1099)

      print*,'TOTAL CPU TIME IS', CPU

      STOP
  125 FORMAT(1X,'LV',6X,' ST',6X,' SV',6X,' SQ',6X,' SC',6X,' SE',6X,'S
     1QA',6X,' UB',6X,' VB',6X,' WB')
  126 FORMAT(1X,'LV',6X,' QC',6X,' QR',6X,' QI',6X,' QS',6X,' QG',6X,'S
     1QC',6X,'SQR',6X,'SQI',6X,'SQS',6X,'SQG')
  130 FORMAT(1X,I2,11(1X,E9.2))
  134 FORMAT(1X,I2,1P14E9.2)

  138 FORMAT(' LV      SMF      SMU      SMD      SMN     CLD      CL
     1DU     CLDD     CLDN     SQF2     SQU2     SQD2     SMN2')
  155 FORMAT(//' THREE-DIMENSIONAL CLOUD ENSEMBLE MODEL'//,
     1  'SECOND ORDER FLUX FORM'///)
 1025 format(1x,'lv',6x,'Q1Z',6x,'Q2Z',6x,'Q1H',6x,'Q2H')
 1099 FORMAT(////5X,'END OF THIS RUN')
12345 FORMAT(////5X,3I6,3X,3F12.2)
66666  format(/5x,'tsfc=',f12.5,4x,'qsfc=',f12.5/)
      END

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ADVECT (X,X1,AA,ISM,IST,ISQ,ISR,ISUV,DLT1,DLT2)
C     ****   COMPUTE ADVECTION OF T, Q, QC, QR, QI, QS AND QG
c
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880,NT2=NT*2)
c
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      common/timestat/ ndt_stat,idq
      COMMON/RTIME/ ISEC1
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1   NRAN,KT1,KT2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(4),AL,CP,RA,CK,CE,EPS,PSFC(5)
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
c
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B2S/ QCS1(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)

      COMMON/B3CR/ QCL2(NX,NY,NZ)

      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/BADV/ DFC(NX,NY,NZ)
      COMMON/BADV1/ VTP(NX,NY,NZ)
      COMMON/BSAT/ XXW(NX,NY,NZ)
      COMMON/BSAT1/ AAA(NX,NY,NZ)
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
c
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSK2(NZ),BSKM(NZ),BSKM4(NZ),
     1   BSI(NX,NY),BSI2(NX,NY),BSIM(NX,NY),BSIM4(NX,NY),
     2   BSJ(NX,NY),BSJ2(NX,NY),BSJM(NX,NY),BSJM4(NX,NY)
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
      COMMON/B4Z/ BSIZ(NZ),BSI2Z(NZ),BSIMZ(NZ),BSIM4Z(NZ),BSJZ(NZ),
     1   BSJ2Z(NZ),BSJMZ(NZ),BSJM4Z(NZ)
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
c
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      COMMON/B7/ SQAAQ(NZ),SQTT(NZ),SQVV(NZ),SQAAT(NZ),SQAK(NZ)

      COMMON/B8/ SMF(NZ),SMU(NZ),SMD(NZ),STF(NZ),STU(NZ),STD(NZ),
     1   SQF(NZ),SQU(NZ),SQD(NZ),SDT(NZ),SDQ(NZ),CLD(NZ),STV(NZ)
      COMMON/B9/ SUT1(NZ),SUC1(NZ),SUN1(NZ),SUU1(NZ),SUD1(NZ),AUB(NZ),
     1   SSU(NZ),SUT2(NZ),SUC2(NZ),SUN2(NZ),SUU2(NZ),SUD2(NZ),AUC1(NZ),
     2   AUN1(NZ),AUU1(NZ),AUD1(NZ),AUC2(NZ),AUN2(NZ),AUU2(NZ),AUD2(NZ),
     3   SVT1(NZ),SVC1(NZ),SVN1(NZ),SVU1(NZ),SVD1(NZ),AVB(NZ),
     4   SSV(NZ),SVT2(NZ),SVC2(NZ),SVN2(NZ),SVU2(NZ),SVD2(NZ),AVC1(NZ),
     5   AVN1(NZ),AVU1(NZ),AVD1(NZ),AVC2(NZ),AVN2(NZ),AVU2(NZ),AVD2(NZ),
     6   STF1(NZ),STU1(NZ),STD1(NZ),SQF1(NZ),SQU1(NZ),SQD1(NZ),
     7   STF2(NZ),STU2(NZ),STD2(NZ),SQF2(NZ),SQU2(NZ),SQD2(NZ),CLDU(NZ),
     8   CLDD(NZ),CLDN(NZ),SMN(NZ),STN1(NZ),STN2(NZ),STN(NZ),SQN1(NZ),
     9   SQN2(NZ),SQN(NZ)
      COMMON/BSFC/ tsfc_1(nx,ny), qsfc_1(nx,ny)
      COMMON/BSTS1/ TUT1(NZ,4,7),TUT2(NZ,4,7),TVT1(NZ,4,7),TVT2(NZ,4,7),
     1 TSTF(NZ,4,7),TSTF1(NZ,4,7),TSTF2(NZ,4,7),TSQF(NZ,4,7),
     2 QQQ(NZ,4,7),TSQF1(NZ,4,7),TSQF2(NZ,4,7),TSQQ(NZ,4,7),
     3 TSQQ1(NZ,4,7)
      COMMON/BSTS2/ snth(nz,4,7),sntv(nz,4,7),sntd(nz,4,7),SNTL(NZ,4,7),
     1 SNQH(NZ,4,7),snqv(nz,4,7),snqd(nz,4,7),SNQL(NZ,4,7),
     2 SNTL1(NZ,4,7),snhh(nz,4,7),snhv(nz,4,7),snhd(nz,4,7),
     3 SNQL1(NZ,4,7)
      COMMON/BSTS3/ QV0(NZ,4,7),TT0(NZ,4,7),SQV0(NZ,4,7),STT0(NZ,4,7),
     1 SGPT(NZ,4,7),SSGPT(NZ,4,7),SNQHD(NZ,4,7),SNQVD(NZ,4,7),
     2 Q1T(NZ,4,7),SNHDH(NZ,4,7),SQHDT(NZ,4,7),SQVDT(NZ,4,7)
      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7),SQTDT(NZ,4,7),SQHL(NZ,4,7)
      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NX,NY),ICS5(NX,NY,4),
     1  IBZ(NX,NY,4)
      COMMON/BSTS20/ OTHERT_ADD(NZ,4,7),OTHERQ_ADD(NZ,4,7)
      common/bsts40/ fcld(nz,4,7)
      COMMON/BTT/ S1(16,NZ),S2(16,NZ),S3(16,NZ),S4(16,NZ),S5(16,NZ),
     1   S6(16,NZ),S7(16,NZ),S8(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     2   S14(16,NZ),S15(16,NZ),S16(16,NZ),SN1(5,NZ),SN2(5,NZ),SN3(5,NZ),
     3   SN4(5,NZ),SN5(5,NZ),SN6(5,NZ),SN7(5,NZ),SN8(5,NZ),SN11(5,NZ),
     4   SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),SN15(5,NZ),SN16(5,NZ)
      common/gbs/ tlsw(nz),qlsw(nz),ttlsw(nz),tqlsw(nz),smfff(nz),
     1   smffff(nz),tcof(nz),textra(nz),qextra(nz)

      COMMON/SURFACE/ SUN_4(NX,NY,4)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON/BBA/ XY1(NX,NY),XY2(NX,NY),XY3(NX,NY),XY4(NX,NY),Y1(NM),
     1   Y2(NM),Y3(NM),Y4(NM),Y5(NM),Y6(NM),Y7(NM),Y8(NM),Y9(NM),
     1   Y10(NM),Y11(NM),Y12(NM),Y13(NM),Y14(NM),Y15(NM),Y16(NM)
      COMMON/RSTAT/ CSTTT(NX,NY),CSTT(NX,NY)
      COMMON/BI/ IT(NX,NY),ICS(NX,NY,4)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
      common/TMP/ y1d(nx,ny,nz),y2d(nx,ny,nz),y3d(nx,ny,nz),
     1            y4d(nx,ny,nz),y03d(nx,ny,1)

      common/bch1/ rby(7)

C   LOCAL VARIABLES

      REAL AA(1)
      REAL X(NX,NY,1)
      REAL X1(NX,NY,1)
      COMMON /tmp1/ XX0(NX,NY,NZ),XX00(NX,NY,NZ),XX000(NX,NY,NZ)
      REAL Y(NX,NY,NZ)
      REAL Y17(NM)
      REAL Y18(NM)
      REAL Y19(NM)
      REAL QQ(NM)

      save

C     ***  STATISTICS FOR CLOUD PROPERTIES   ***********************

      BARFCT=2.*100.*RHO1(2)    !CONVERT FROM MKS TO CGS
      D22T=D2T
      assss=1.
      awks=1.
      d2t=d2t
      IF (IJKADV .EQ. 1) THEN
         D2T=DT
         assss=0.
         awks=1.
      ENDIF

      a0000=1.
      if (imlifting .eq. 0) a0000=0.

      a_irf=1.
      if (irf .eq. 0) a_irf=0.

      a0 =ndt_stat*0.5*rijl2
      a20=ndt_stat*rijl2
      a300=d2t*rijl2
C
      IF(ISM+IST+ISQ+IDQ .NE. 0) THEN
         if (idq .eq. 0) then

            IF(ISM.NE.0) THEN
c$doacross local(k,km,yy7,yy8,yy15,j,jp,i,ip,yy1,yy5,yy6,yy9,yy10,yy11,
c$&             yy12,yy13,yy14,yy17,yy18,qqqq,ww,rq,iww,jww,kc,mt,ibz1)
               DO K=2,KLES
                  KM=K-1
                  YY7=.5*(UB1(K)+UB1(KM))
                  YY8=.5*(VB1(K)+VB1(KM))
                  DO J=2,JLES
                     JP=J+1
                     DO I=2,ILES
                        IP=I+1
                        YY1=RHO1(K)*(WW1(I,J,K)+y0(I,j)*WB(K))
                        YY5=.25*(UU1(I,J,K)+UU1(IP,J,K)+UU1(I,J,KM)
     1                     +UU1(IP,J,KM))
                        YY6=.25*(VV1(I,J,K)+VV1(I,JP,K)+VV1(I,J,KM)
     1                     +VV1(I,JP,KM))
                        YY9=YY5-YY7
                        YY10=YY6-YY8
                        YY11=YY1*YY9
                        YY12=YY1*YY10
                        YY13=YY1*YY5
                        YY14=YY1*YY6
                        yy15=.01*(ww1(i,j,k)+wb(k)*y0(i,j))
                        YY17=QCL(I,J,K)+QCI(I,J,K)
                        YY18=QCL(I,J,KM)+QCI(I,J,KM)
                        QQQQ=YY17+YY18+QRN(I,J,K)+QRN(I,J,KM)+QCS(I,J,K)
     1                      +QCG(I,J,K)+QCS(I,J,KM)+QCG(I,J,KM)
                        WW=2.*A0*YY1
                        rq=.01*(ww1(i,j,k)+wb(k))
                        iww=max(0,min(int(rq+0.5), 15) + 1)
                        jww=abs(min(0,max(int(rq-0.5), -5)))
                        SMF(K)=SMF(K)+WW
                        SUT1(K)=SUT1(K)+YY11
                        SVT1(K)=SVT1(K)+YY12
                        SUT2(K)=SUT2(K)+YY13
                        SVT2(K)=SVT2(K)+YY14
                        IF(QQQQ.GE.2.E-5) THEN
                           CLD(K)=CLD(K)+1.
                           SUC1(K)=SUC1(K)+YY11
                           SUC2(K)=SUC2(K)+YY13
                           AUC1(K)=AUC1(K)+YY5
                           AUC2(K)=AUC2(K)+YY9
                           SVC1(K)=SVC1(K)+YY12
                           SVC2(K)=SVC2(K)+YY14
                           AVC1(K)=AVC1(K)+YY6
                           AVC2(K)=AVC2(K)+YY10
                           IF(WW.GE.0.) THEN
                              CLDU(K)=CLDU(K)+1.
                              SMU(K)=SMU(K)+WW
                              SUU1(K)=SUU1(K)+YY11
                              SUU2(K)=SUU2(K)+YY13
                              AUU1(K)=AUU1(K)+YY5
                              AUU2(K)=AUU2(K)+YY9
                              SVU1(K)=SVU1(K)+YY12
                              SVU2(K)=SVU2(K)+YY14
                              AVU1(K)=AVU1(K)+YY6
                              AVU2(K)=AVU2(K)+YY10
                              S1(IWW,K)=S1(IWW,K)+WW
                              S2(IWW,K)=S2(IWW,K)+1.
                              S11(IWW,K)=S11(IWW,K)+YY9
                              S12(IWW,K)=S12(IWW,K)+YY10
                              S13(IWW,K)=S13(IWW,K)+YY5
                              S14(IWW,K)=S14(IWW,K)+YY6
                              S15(IWW,K)=S15(IWW,K)+YY11
                              S16(IWW,K)=S16(IWW,K)+YY12
                           ELSE
                              CLDD(K)=CLDD(K)+1.
                              SMD(K)=SMD(K)+WW
                              SUD1(K)=SUD1(K)+YY11
                              SUD2(K)=SUD2(K)+YY13
                              AUD1(K)=AUD1(K)+YY5
                              AUD2(K)=AUD2(K)+YY9
                              SVD1(K)=SVD1(K)+YY12
                              SVD2(K)=SVD2(K)+YY14
                              AVD1(K)=AVD1(K)+YY6
                              AVD2(K)=AVD2(K)+YY10
                              IF(JWW.NE.0) THEN
                                 SN1(JWW,K)=SN1(JWW,K)+WW
                                 SN2(JWW,K)=SN2(JWW,K)+1.
                                 SN11(JWW,K)=SN11(JWW,K)+YY9
                                 SN12(JWW,K)=SN12(JWW,K)+YY10
                                 SN13(JWW,K)=SN13(JWW,K)+YY5
                                 SN14(JWW,K)=SN14(JWW,K)+YY6
                                 SN15(JWW,K)=SN15(JWW,K)+YY11
                                 SN16(JWW,K)=SN16(JWW,K)+YY12
                              ENDIF
                           ENDIF
                        ELSE
                           CLDN(K)=CLDN(K)+1.
                           SMN(K)=SMN(K)+WW
                           SUN1(K)=SUN1(K)+YY11
                           SUN2(K)=SUN2(K)+YY13
                           AUN1(K)=AUN1(K)+YY5
                           AUN2(K)=AUN2(K)+YY9
                           SVN1(K)=SVN1(K)+YY12
                           SVN2(K)=SVN2(K)+YY14
                           AVN1(K)=AVN1(K)+YY6
                           AVN2(K)=AVN2(K)+YY10
                        ENDIF
                        DO KC=1,7
                           DO MT=1,4
                     IBZ1=ICS5(I,J,MT)
                     IF(KC.GT.4 .and. YY15.GT.RBY(KC))IBZ1=0
                     IF(KC.LT.4 .and. YY15.LT.RBY(KC))IBZ1=0
                     IF(MT.EQ.1) IBZ1=1
                           tut1(k,mt,kc)=tut1(k,mt,kc)+yy11*ibz1
                           tut2(k,mt,kc)=tut2(k,mt,kc)+yy13*ibz1
                           tvt1(k,mt,kc)=tvt1(k,mt,kc)+yy12*ibz1
                           tvt2(k,mt,kc)=tvt2(k,mt,kc)+yy14*ibz1
                       if(qqqq.ge.2.e-5)fcld(k,mt,kc)=fcld(k,mt,kc)+ibz1
                           ENDDO!MT
                        ENDDO	!KC
                     ENDDO	!I
                  ENDDO		!J
               ENDDO		!K
            ENDIF

            IF(IST.NE.0) THEN
c$doacross local(k,km,yy15,aakm,a0k,j,i,yy1,yy3,yy4,yy17,yy18,yyy16,
c$&              yy16,ww,rq,iww,jww,qqqq)
               DO K=2,KLES
                  KM=K-1
                  YY15=FLOAT(ISTATMIN)*RIJL2*RDZ*AM1(K)*RRHO(K)
                  aakm=aa(k)-aa(km)
                  A0K=awks*am1(k)*bsk2(k)
                  DO J=2,JLES
                     DO I=2,ILES
                        YY1=RHO1(K)*(WW1(I,J,K)+y0(I,j)*WB(K))
                        YY3=.5*(X(I,J,K)+X(I,J,KM))
                        YY4=YY1*YY3
                        YY17=QCL(I,J,K)+QCI(I,J,K)
                        YY18=QCL(I,J,KM)+QCI(I,J,KM)
                        IF (YY17 .LT. 1.E-5 .OR. YY18 .LT. 1.E-5) THEN
                         yyy16=((ak(i,j,k)*(1.+xxw(i,j,k))+ak(i,j,km)
     2                   *(1.+xxw(i,j,km)))*am1(k)*(x1(i,j,k)-x1(i,j,km)
     3                   +aakm)+a0k*(x1(i,j,k)-x1(i,j,km)))
                         YY16=rho1(k)*(2.*assss*(ww1(i,j,k)+wb(k)
     1                       *y0(i,j))*yy3-yyy16*rdz)*a0
                        else
                           YY16=rho1(k)*(2.*assss*(ww1(i,j,k)+wb(k)
     1                      *y0(i,j))*yy3+dfc(i,j,k)-a0k*(x1(i,j,k)
     1                      -x1(i,j,km))*rdz)*a0
                        ENDIF
                        WW=2.*A0*YY1
                        rq=.01*(ww1(i,j,k)+wb(k))
                        iww=max(0,min(int(rq+0.5), 15) + 1)
                        jww=abs(min(0,max(int(rq-0.5), -5)))
                        QQQQ=YY17+YY18+QRN(I,J,K)+QRN(I,J,KM)+QCS(I,J,K)
     1                      +QCG(I,J,K)+QCS(I,J,KM)+QCG(I,J,KM)

                        STF1(K)=STF1(K)+YY4
                        STF2(K)=STF2(K)+YY3
                        STF(K)=STF(K)+YY16
                        SDT(K)=SDT(K)+YY1*AAKM*YY15

                        IF(QQQQ.GE.2.E-5) THEN
                           IF(WW.GE.0.) THEN
                              STU1(K)=STU1(K)+YY4
                              STU2(K)=STU2(K)+YY3
                              STU(K)=STU(K)+YY16
                              S3(IWW,K)=S3(IWW,K)+YY4
                              S4(IWW,K)=S4(IWW,K)+YY3
                              S5(IWW,K)=S5(IWW,K)+YY16
                           ELSE 
                              STD1(K)=STD1(K)+YY4
                              STD2(K)=STD2(K)+YY3
                              STD(K)=STD(K)+YY16

                              IF(JWW.NE.0) THEN
                                 SN3(JWW,K)=SN3(JWW,K)+YY4
                                 SN4(JWW,K)=SN4(JWW,K)+YY3
                                 SN5(JWW,K)=SN5(JWW,K)+YY16
                              ENDIF
                           ENDIF
                        ELSE
                           STN1(K)=STN1(K)+YY4
                           STN2(K)=STN2(K)+YY3
                           STN(K)=STN(K)+YY16
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ELSEIF(ISQ.NE.0) THEN
c$doacross local(k,km,yy15,aakm,a0k,j,i,yy1,yy3,yy4,yy17,yy18,yyy16,
c$&              yy16,ww,rq,iww,jww,qqqq)
               DO K=2,KLES
                  KM=K-1
                  YY15=FLOAT(ISTATMIN)*RIJL2*RDZ*AM1(K)*RRHO(K)
                  KM=K-1
                  aakm=aa(k)-aa(km)
                  A0K=awks*am1(k)*bsk2(k)
                  DO J=2,JLES
                     DO I=2,ILES
                        YY1=RHO1(K)*(WW1(I,J,K)+y0(I,j)*WB(K))
                        YY3=.5*(X(I,J,K)+X(I,J,KM))
                        YY4=YY1*YY3
                        YY17=QCL(I,J,K)+QCI(I,J,K)
                        YY18=QCL(I,J,KM)+QCI(I,J,KM)
                        IF (YY17 .LT. 1.E-5 .OR. YY18 .LT. 1.E-5) THEN
                         yyy16=((ak(i,j,k)*(1.+xxw(i,j,k))+ak(i,j,km)
     2                   *(1.+xxw(i,j,km)))*am1(k)*(x1(i,j,k)-x1(i,j,km)
     3                   +aakm)+a0k*(x1(i,j,k)-x1(i,j,km)))
                         YY16=rho1(k)*(2.*assss*(ww1(i,j,k)+wb(k)
     1                       *y0(i,j))*yy3-yyy16*rdz)*a0
                        else
                          YY16=rho1(k)*(2.*assss*(ww1(i,j,k)+wb(k)
     1                      *y0(i,j))*yy3+dfc(i,j,k)-a0k*(x1(i,j,k)
     1                      -x1(i,j,km))*rdz)*a0
                        ENDIF
                        WW=2.*A0*YY1
                        rq=.01*(ww1(i,j,k)+wb(k))
                        iww=max(0,min(int(rq+0.5), 15) + 1)
                        jww=abs(min(0,max(int(rq-0.5), -5)))
                        QQQQ=YY17+YY18+QRN(I,J,K)+QRN(I,J,KM)+QCS(I,J,K)
     1                      +QCG(I,J,K)+QCS(I,J,KM)+QCG(I,J,KM)

                        SQF1(K)=SQF1(K)+YY4
                        SQF2(K)=SQF2(K)+YY3
                        SQF(K)=SQF(K)+YY16
                        SDQ(K)=SDQ(K)+YY1*AAKM*YY15

                        IF(QQQQ.GE.2.E-5) THEN
                           IF(WW.GE.0.) THEN
                              SQU1(K)=SQU1(K)+YY4
                              SQU2(K)=SQU2(K)+YY3
                              SQU(K)=SQU(K)+YY16
                              S6(IWW,K)=S6(IWW,K)+YY4
                              S7(IWW,K)=S7(IWW,K)+YY3
                              S8(IWW,K)=S8(IWW,K)+YY16
                           ELSE 
                              SQD1(K)=SQD1(K)+YY4
                              SQD2(K)=SQD2(K)+YY3
                              SQD(K)=SQD(K)+YY16
                              IF(JWW.NE.0) THEN
                                 SN6(JWW,K)=SN6(JWW,K)+YY4
                                 SN7(JWW,K)=SN7(JWW,K)+YY3
                                 SN8(JWW,K)=SN8(JWW,K)+YY16
                              ENDIF
                           ENDIF
                        ELSE
                           SQN1(K)=SQN1(K)+YY4
                           SQN2(K)=SQN2(K)+YY3
                           SQN(K)=SQN(K)+YY16
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

            DO K=2,KLES
               KP=K+1
               KM=K-1
               WBKP=WB(KP)+WB(K)
               YY15=FLOAT(ISTATMIN)*RIJL2*RDZ*AM1(K)*RRHO(K)
               Y7(K)=.5*(UB1(K)+UB1(KM))
               Y8(K)=.5*(VB1(K)+VB1(KM))
               a000c=a0000*a20*am(k)
               A0K=awks*am1(k)*bsk2(k)
               a_irf1=a_irf*a20*RFA(K)
               aakm=aa(k)-aa(km)
               DO J=2,JLES
                  JP=J+1
                  DO 600 I=2,ILES
                     IP=I+1
                     YY1=RHO1(K)*(WW1(I,J,K)+y0(I,j)*WB(K))
                     YY3=.5*(X(I,J,K)+X(I,J,KM))
                     YY5=.25*(UU1(I,J,K)+UU1(IP,J,K)+UU1(I,J,KM)
     1                  +UU1(IP,J,KM))
                     YY6=.25*(VV1(I,J,K)+VV1(I,JP,K)+VV1(I,J,KM)
     1                  +VV1(I,JP,KM))
                     YY9=YY5-Y7(K)
                     YY10=YY6-Y8(K)
                     Y11(i)=YY1*YY9
                     Y12(i)=YY1*YY10
                     Y13(i)=YY1*YY5
                     Y14(i)=YY1*YY6
                     YY17=QCL(I,J,K)+QCI(I,J,K)
                     YY18=QCL(I,J,KM)+QCI(I,J,KM)
                     QQ(i)=YY17+YY18+QRN(I,J,K)+QRN(I,J,KM)+QCS(I,J,K)
     1                    +QCG(I,J,K)+QCS(I,J,KM)+QCG(I,J,KM)
                  yy20=0.0 
                  if (k.eq.2) then
                     IF (ISmg.EQ.1) YY20=tsfc_1(i,j)
                     IF (ISmg.EQ.2) YY20=qsfc_1(i,j)
                     IF (ITOGA .EQ. 1) THEN
                       if (ismg .eq. 1) YY20=100.*rho1(2)*swt(i,j)/pi(2)
                       if (ismg .eq. 2) YY20=100.*rho1(2)*swq(i,j)
                     ENDIF
                  endif
                  Y19(I)=-a000c*y0(i,j)*WBKP*(X(I,J,KP)
     1                  -X(I,J,KM))*RD4Z-a_irf1*X1(I,J,K)
               if(yy17.lt.1.e-5 .or. yy18.lt.1.e-5) then
                       yyy16=((ak(i,j,k)*(1.+xxw(i,j,k))+ak(i,j,km)
     2                   *(1.+xxw(i,j,km)))*am1(k)*(x1(i,j,k)-x1(i,j,km)
     3                   +aakm)+a0k*(x1(i,j,k)-x1(i,j,km)))
                  y16(i)=rho1(k)*(2.*assss*ww1(i,j,k)*yy3
     1                 -yyy16*rdz)*a0+yy20*a20
               else
                  y16(i)=rho1(k)*(2.*assss*ww1(i,j,k)*yy3+dfc(i,j,k)
     1                  -a0k*(x1(i,j,k)-x1(i,j,km))*rdz)*a0+yy20*a20
               endif
               y15(i)=.01*(ww1(i,j,k)+wb(k)*y0(i,j))
  600      continue
            DO 900 KC=1,7
               DO 98 MT=1,4
                  DO I=2,ILES
                     IBZ(I,J,MT)=ICS5(I,J,MT)
                     IF(KC.GT.4 .and. Y15(I).GT.RBY(KC)) IBZ(I,J,MT)=0
                     IF(KC.LT.4 .and. Y15(I).LT.RBY(KC)) IBZ(I,J,MT)=0
                     IBZ(I,J,1)=1
                  enddo
                  if(ist.eq.1) then
                     A5=0.
                     A6=0.
                     A7=0.
                     A55=0.
                    do I=2,ILES
                        yy21= (x(i,j,k)+ta(k)-x1(i,j,k)-ta1(k))/dt
                        yy22=(am1(k)*fd(k)*ww1(i,j,k)
     1                      +am1(kp)*fd(kp)*ww1(i,j,kp))*rdz*rrho(k)
                        a5=a5+y16(i)*ibz(i,j,mt)
                        a6=a6+yy21*a20*ibz(i,j,mt)
                        a7=a7+yy22*a20*ibz(i,j,mt)
                        A55=A55+Y19(I)*ibz(i,j,mt)
                     enddo
                     tstf(k,mt,kc)=tstf(k,mt,kc)+a5
                     tstf2(k,mt,kc)=tstf2(k,mt,kc)+a6
                     tstf1(k,mt,kc)=tstf1(k,mt,kc)+a7
                     othert_add(k,mt,kc)=othert_add(k,mt,kc)+a55
                  endif
                  if(isq.eq.1) then
                     A8=0.
                     A9=0.
                     A10=0.
                     A66=0.
                     DO I=2,ILES
                        yy21=(x(i,j,k)+qa(k)-x1(i,j,k)-qa1(k))/dt
                        yy22=(am1(k)*fe(k)*ww1(i,j,k)
     1                      +am1(kp)*fe(kp)*ww1(i,j,kp))*rdz*rrho(k)
                        a8=a8+y16(i)*ibz(i,j,mt)
                        a9=a9+yy21*a20*ibz(i,j,mt)
                        a10=a10+yy22*a20*ibz(i,j,mt)
                        A66=A66+Y19(I)*ibz(i,j,mt)
                     enddo
                     tsqf(k,mt,kc)=tsqf(k,mt,kc)+a8
                     tsqf2(k,mt,kc)=tsqf2(k,mt,kc)+a9
                     tsqf1(k,mt,kc)=tsqf1(k,mt,kc)+a10
                     otherq_add(k,mt,kc)=otherq_add(k,mt,kc)+a66
                  endif
   98          continue
  900       continue
           ENDDO
           ENDDO
      endif

C     ****** STATISTIC FOR U- AND V-MOMENTUM

      IF(ISUV.EQ.1) THEN
         DO K=2,KLES
            AUB(K)=AUB(K)+Y7(K)
            AVB(K)=AVB(K)+Y8(K)
         ENDDO
      ENDIF

c     ****** statistic for w(qc+qr+qi+qs+qg)
      if (idq.eq.1) then
         do 7001 k=2,kles
            km=k-1
            A0K=AWKS*AM1(K)*BSK2(K)
            do 7000 j=2,jles
               do 71 i=2,iles
                  Y3(I)=RHO1(K)*(assss*ww1(i,j,k)*(X(I,J,K)+X(I,J,KM))
     1                 -(AM1(K)*(X1(I,J,K)-X1(I,J,KM))
     2                 *(AK(I,J,K)*(1.+XXW(I,J,K))+AK(I,J,KM)*(1.
     3                +XXW(I,J,KM)))+A0K*(X1(I,J,K)-X1(I,J,KM)))*RDZ)*a0
                  IF(ISR.NE.1) THEN
                     IF (ISMG .EQ. 4) THEN
                        Y17(I)=DLT1*QCL(I,J,K)+DLT2*QCI(I,J,K)
                        Y18(I)=DLT1*QCL(I,J,KM)+DLT2*QCI(I,J,KM)
                     ELSE
                        Y17(I)=DLT1*QCL2(I,J,K)+DLT2*QCI(I,J,K)
                        Y18(I)=DLT1*QCL2(I,J,KM)+DLT2*QCI(I,J,KM)
                     ENDIF
                     if(y17(i).GE.1.e-6 .AND. y18(i).GE.1.e-6) THEN
                        y3(i)=rho1(k)*(assss*ww1(i,j,k)*(x(i,j,k)
     1                          -x(i,j,km))+dfc(i,j,k)-a0k
     2                          *(x1(i,j,k)-x1(i,j,km))*rdz)*a0
                     endif
                  ENDIF
   71          continue
               do 700 kc=1,7
                  do mt=1,4
                     do i=2,iles
                        y2(i)=.01*ww1(i,j,k)
                        ibz(i,j,mt)=0
                        ibz(i,j,1)=1
                        if(ics5(i,j,mt).eq.1) ibz(i,j,mt)=1
                        if(kc.gt.4) then
                           if (y2(i).gt.rby(kc)) ibz(i,j,mt)=0
                        elseif(kc.lt.4)then
                           if (y2(i).lt.rby(kc)) ibz(i,j,mt)=0
                        endif
                   enddo
                 enddo
                  do mt=1,4
                     a11=0.
                     a12=0.
                     do i=2,iles
                        yy21=(x(i,j,k)-x1(i,j,k))/dt
                        a11=a11+yy21*a20*ibz(i,j,mt)
                        a12=a12+y3(i)*ibz(i,j,mt)
                     enddo
                     tsqq(k,mt,kc)=tsqq(k,mt,kc)+a11
                     tsqq1(k,mt,kc)=tsqq1(k,mt,kc)+a12
                  enddo
  700          continue
 7000    continue
 7001    continue

c     ******   vertical advection terms   **************************
      endif
      endif
c
      if (ismg .eq. 1) then
cc$doacross local(k,j,i,aextra)
         do k=2,kles
            do j=2,jles
               do i=2,iles
                  aextra=am(k)*y0(i,j)*(wb(k+1)+wb(k))
     1                  *(x(i,j,k+1)-x(i,j,k-1))*rd4z*a300
                   textra(k)=textra(k)-aextra
               enddo
            enddo
         enddo
      elseif (ismg .eq. 2) then
cc$doacross local(k,j,i,aextra)
         do k=2,kles
            do j=2,jles
               do i=2,iles
                  aextra=am(k)*y0(i,j)*(wb(k+1)+wb(k))
     1                  *(x(i,j,k+1)-x(i,j,k-1))*rd4z*a300
                   qextra(k)=qextra(k)-aextra
               enddo
            enddo
         enddo
      endif

C-----put surface layer heat or moisture fluxs in
c
      if (ijkadv .eq. 0) then
c
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
        XY3(I,J)=0.0
        XY1(I,J)=RHO1(2)*WW1(I,J,2)*(X(I,J,2)+X(I,J,1))
      ENDDO
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     *****     TOGA-COARE SENSIBLE AND LATENT HEAT FLUX    ************
      IF (ITOGA .EQ. 1) THEN
c$doacross local(j,i)
        DO J=2,JLES
        DO I=2,ILES
            if (ismg .eq. 1) then
              xy3(i,j)=swt(i,j)*barfct/pi(2)
              sun_4(i,j,3)=100.*swt(i,j)/pi(2)
            endif
            if (ismg .eq. 2) then
              xy3(i,j)=swq(i,j)*barfct
              sun_4(i,j,4)=100.*swq(i,j)
            endif
        ENDDO
        ENDDO
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 1012 K=3,KMAX
       Y1(K)=AM1(K)*BSK2(K)
 1012 CONTINUE
      DO 10001 K=3,KMAX
        KM=K-1
        A1=Y1(K)*RDZ
        A0K=AM1(K)*BSK2(K)
        A1K=A0K*RDZ
        A2K=AM(KM)*RRHO(KM)*RD2Z

c$doacross local(j,i)
      DO J=2,JLES
         DO I=2,ILES
            XY2(I,J)=RHO1(K)*WW1(I,J,K)*(X(I,J,K)+X(I,J,KM))
         ENDDO
      ENDDO

c$doacross local(j,i)
      DO J=2,JLES
        DO I=2,ILES
         XY4(I,J)=RHO1(K)*(-(AM1(K)*(X1(I,J,K)-X1(I,J,KM)+AA(K)-AA(KM))
     1      *(AK(I,J,K)*(1.+XXW(I,J,K))+AK(I,J,KM)*(1.+XXW(I,J,KM)))
     2       +Y1(K)*(X1(I,J,K)-X1(I,J,KM)))*RDZ)
        ENDDO
      ENDDO

      DO 10000 J=2,JLES
       IF(ISR.EQ.0) THEN
        DO 1400 I=2,ILES
          IF (ISMG .EQ. 4) THEN
           Y3(I)=DLT1*QCL2(I,J,K)+DLT2*QCI(I,J,K)
           Y4(I)=DLT1*QCL2(I,J,KM)+DLT2*QCI(I,J,KM)
          ELSE
           Y3(I)=DLT1*QCL(I,J,K)+DLT2*QCI(I,J,K)
           Y4(I)=DLT1*QCL(I,J,KM)+DLT2*QCI(I,J,KM)
          ENDIF
         if(y3(i) .lt. 1.e-5 .or. y4(i) .lt. 1.e-5) go to 1400
           xy4(i,j)=rho1(k)*(dfc(i,j,k)-a1k*(x1(i,j,k)-x1(i,j,km)))
 1400   CONTINUE
       ENDIF

       if (k.eq.kmax) then  !shn
          do i=2,iles
             xy2(i,j)=0.
             xy4(i,j)=0.
          enddo
       endif
       DO 1500 I=2,ILES
       Y(I,J,KM)=AM(KM)*(XY1(I,J)-XY2(I,J)+XY3(I,J)-XY4(I,J))
     1           *RRHO(KM)*RD2Z
c
ctao
         XY1(I,J)=XY2(I,J)
         XY3(I,J)=XY4(I,J)
 1500   continue
10000  CONTINUE
10001  CONTINUE

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-18-97)

          r2dt8=0.008/d2t

          y1(2)=0.
          y1(kmax)=0.

cc$doacross local(j,i)
          do j=2,jles
             do i=2,iles
                yyy13=x1(i,j,3)-x1(i,j,2)
                y(i,j,2)=y(i,j,2)-r2dt8*(y1(2)-yyy13)
             enddo
          enddo

c$doacross local(k,km,j,i)
          do k=3,kles
             km=k-1
             do j=2,jles
                do i=2,iles
                   y1(k)=x1(i,j,k)-x1(i,j,km)
                   y(i,j,k)=y(i,j,k)-r2dt8*(y1(k)-y1(k+1))
                enddo
             enddo
            enddo


cccccc  tao (11-18-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     ******   HORIZONTAL ADVECTION   ******


      DO 30000 K=2,KLES
       KP=K+1
       KM=K-1
       amkk=AM(K)*RD4Z

       IF (IADVH .EQ. 2) THEN
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,IMAX
        XY1(I,J)=UU1(I,J,K)*(X(I,J,K)+X(I-1,J,K))*RD2X
       ENDDO
       ENDDO

c$doacross local(j,jm,i)
       DO J=2,JMAX
       JM=J-1
       DO I=2,ILES
        XY3(I,J)=VV1(I,J,K)*(X(I,J,K)+X(I,JM,K))*RD2Y
       ENDDO
       ENDDO

c$doacross local(j,jp,i)
       DO J=2,JLES
       JP=J+1
       DO I=2,ILES
        Y(I,J,K)=Y(I,J,K)+(XY1(I,J)-XY1(I+1,J))
     1                   +(XY3(I,J)-XY3(I,JP))
       ENDDO
       ENDDO
      ELSE
CC    ******   4-TH ORDER  HORIZONTAL ADVECTION TERMS   ****************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c$doacross local(j,i,ipp,ip)
       DO J=2,JLES
        do i=2,iles
           ipp=i+2
           if (i.eq.iles) ipp=3
           ip=i1
         xy2(i,j)=uu1(ip,j,k)*(d58x*(x(i,j,k)+x(ip,j,k))
     1                     -d16x*(x(i-1,j,k)+x(ipp,j,k)))
       enddo
       enddo

c$doacross local(j,i)
       DO J=2,JLES
        do i=2,imax
         xy3(i,j)=-d48x*uu1(i,j,k)*(x(i,j,k)+x(i-1,j,k))
        enddo
        enddo

c$doacross local(j,i,im,ipp)
       DO 3040 J=2,JLES
        do 306 i=2,iles
           im=i-1
             if(i.eq.2) im=iles !!sh
           ipp=i+2
             if (i.eq.iles) ipp=3
         y(i,j,k)=y(i,j,k)+(xy2(im,j)-xy2(i,j)+xy3(im,j)-xy3(ipp,j))
  306    continue
 3040  CONTINUE
CC   ****************

c$doacross local(j,jm,i)
      do j=2,jmax
         jm=j-1
         DO I=2,ILES
            xy3(i,j)=-d48y*vv1(i,j,k)*(x(i,j,k)+x(i,jm,k))
         enddo
      enddo

c$doacross local(j,jpp,jp,jm,i)
       do j=2,jles
          jpp=j+2
          if(j.eq.jles) jpp=3
          jp=j+1
          jm=j-1
          DO I=2,ILES
          xy1(i,j)=vv1(i,jp,k)*(d58y*(x(i,j,k)+x(i,jp,k))
     1                          -d16y*(x(i,jm,k)+x(i,jpp,k)))
        enddo
        enddo

c$doacross local(j,jm,jpp,i)
       do j=2,jles
         jm=j-1
         if(j.eq.2) jm=jles
         jpp=j+2
         if(j.eq.jles) jpp=3
         DO I=2,ILES
         y(i,j,k)=y(i,j,k)+(xy1(i,jm)-xy1(i,j)+xy3(i,jm)-xy3(i,jpp))
        enddo
        enddo
C     ******************************************************************
       ENDIF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$doacross local(j,i)
       DO J=2,JLES
        DO I=2,IMAX
         XY1(I,J)=-R2DX2*(AK(I,J,K)+AK(I-1,J,K))*(X1(I,J,K)-X1(I-1,J,K))
        ENDDO
        ENDDO

c$doacross local(j,jm,i)
       DO J=2,JMAX
       JM=J-1
        DO I=2,ILES
         XY2(I,J)=-R2DY2*
     1            (AK(I,J,K)+AK(I,JM,K))*(X1(I,J,K)-X1(I,JM,K))
        ENDDO
        ENDDO

c$doacross local(j,jp,i)
       DO J=2,JLES
       JP=J+1
        DO I=2,ILES
         Y(I,J,K)=Y(I,J,K)+(XY1(I,J)-XY1(I+1,J))+(XY2(I,J)-XY2(I,JP))

        ENDDO
        ENDDO
c

c$doacross local(j,i,ipp)
       DO J=2,JLES
        do i=2,iles
           ipp=i+2
             if (i.eq.iles) ipp=3
             xy2(i,j)=a4h(k)*(x1(ipp,j,k)-3.*(x1(i+1,j,k)-x1(i,j,k))
     1                   -x1(i-1,j,k))
         enddo
         enddo

c$doacross local(j,i,im)
       DO 2400 J=2,JLES
        do  i=2,iles
           im=i-1
           if(i.eq.2) im=iles 
         y(i,j,k)=y(i,j,k)+(xy2(im,j)-xy2(i,j))
        enddo
 2400  CONTINUE

c$doacross local(j,jp,jm,jpp,i)
       do j=2,jles
       jp=j+1
       jm=j-1
       jpp=j+2
       if(j.eq.jles) jpp=3
       DO I=2,ILES
         xy1(i,j)=a4h(k)*(x1(i,jpp,k)-3.*(x1(i,jp,k)-x1(i,j,k))
     1                   -x1(i,jm,k))
         enddo
         enddo

cdaniel
c$doacross local(j,i,jm)
       do j=2,jles
       DO I=2,ILES
           jm=j-1
           if(j.eq.2) jm=jles !!sh
         y(i,j,k)=y(i,j,k)+(xy1(i,jm)-xy1(i,j))
        enddo
        enddo

c$doacross local(j,i)
       DO 2081 J=2,JLES
       DO 2080 I=2,ILES
           Y(I,J,K)=Y(I,J,K)-a0000*AM(K)*Y0(I,J)*(WB(KP)+WB(K))*
     1           (X(I,J,KP)-X(I,J,KM))*RD4Z-a_irf*RFA(K)*X1(I,J,K)
 2080  CONTINUE 
 2081  CONTINUE 
30000 CONTINUE
C

      DO K=2,KLES
c$doacross local(j,i)
         DO J=2,JLES
            DO I=2,ILES
               XY1(I,J)=X(I,J,K)+EPS*(-2.*X(I,J,K)+X1(I,J,K))
               X(I,J,K)=X1(I,J,K)+Y(I,J,K)*D2T
               X1(I,J,K)=XY1(I,J)
            ENDDO   
         ENDDO
      ENDDO
       
CC
c     ****************************************************************

      ELSE  !!sh (if (ijkadv .ne. 0))!

c     ****************************************************************
c$doacross local(j,i)
      DO J=2,JLES
         DO I=2,ILES
            XY3(I,J)=0.0
         ENDDO
      ENDDO

C     *****     TOGA-COARE SENSIBLE AND LATENT HEAT FLUX    ************
      IF (ITOGA .EQ. 1) THEN

            if (ismg .eq. 1) then
c$doacross local(j,i)
        DO J=2,JLES
        DO I=2,ILES
              xy3(i,j)=swt(i,j)*barfct/pi(2)
              sun_4(i,j,3)=100.*swt(i,j)/pi(2)
        ENDDO
        ENDDO
            endif
            if (ismg .eq. 2) then
c$doacross local(j,i)
        DO J=2,JLES
        DO I=2,ILES
              xy3(i,j)=swq(i,j)*barfct
              sun_4(i,j,4)=100.*swq(i,j)
        ENDDO
        ENDDO
            endif
      ENDIF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      BWK=awks


c$doacross local(k,rho1k,am1k,km,aak_aakm,a0k,j,i)
      DO K=3,KMAX
         RHO1K=RHO1(K)
         AM1K=AM1(K)
         KM=K-1
c         AAK=AA(K)
c         AAKM=AA(K-1)
         aak_aakm=AA(K)-AA(K-1)
         A0K=awks*AM1(K)*BSK2(K)
         DO J=2,JLES
            DO I=2,ILES
              Y4D(I,J,K)=RHO1K*(-(AM1K*(X1(I,J,K)-X1(I,J,KM)+aak_aakm)
     1          *(AK(I,J,K)*(1.+XXW(I,J,K))+AK(I,J,KM)*(1.+XXW(I,J,KM)))
     2          +A0K*(X1(I,J,K)-X1(I,J,KM)))*RDZ)
            ENDDO
         ENDDO
      ENDDO

      IF(ISR.EQ.0) THEN
         DO K=3,KMAX
         KM=K-1
         A1K=(awks*AM1(K)*BSK2(K))*RDZ
cc         A2K=AM(KM)*RRHO(KM)*RD2Z
         DO J=2,JLES
            DO I=2,ILES
               IF (ISMG .EQ. 4) THEN
                  Y3(I)=DLT1*QCL2(I,J,K)+DLT2*QCI(I,J,K)
                  Y4(I)=DLT1*QCL2(I,J,KM)+DLT2*QCI(I,J,KM)
               ELSE
                  Y3(I)=DLT1*QCL(I,J,K)+DLT2*QCI(I,J,K)
                  Y4(I)=DLT1*QCL(I,J,KM)+DLT2*QCI(I,J,KM)
               ENDIF
               if(y3(i) .LT. 1.e-5 .OR. y4(i) .LT. 1.e-5) GO TO 140
              y4d(i,j,k)=RHO1(K)*(dfc(i,j,k)-a1k*(x1(i,j,k)-x1(i,j,km)))
  140          CONTINUE
            ENDDO
         ENDDO
         ENDDO
      ENDIF
c
c$doacross local(j,i)
      do j=2,jles
      do i=2,iles
         y4d(i,j,kmax)=0.
      enddo
      enddo
c
      DO K=3,KMAX
         KM=K-1
cc         A1K=(awks*AM1(K)*BSK2(K))*RDZ
         A2K=AM(KM)*RRHO(KM)*RD2Z
c$doacross local(j,i)
         DO J=2,JLES
            DO I=2,ILES
               Y(I,J,KM)=A2K*(XY3(I,J)-Y4D(I,J,K))
               XY3(I,J)=Y4D(I,J,K)
            ENDDO
         ENDDO
      ENDDO

c$doacross local(k,j,i)
      DO K=2,KLES
         DO J=2,JLES
            DO I=2,ILES
               XX000(I,J,K)=Y(I,J,K)
            ENDDO
         ENDDO
      ENDDO

C     ******   HORIZONTAL ADVECTION   *********************************


c$doacross local(k,j,i)
      DO K=2,KLES
         DO J=2,JLES
            DO I=2,IMAX
               Y1D(I,J,K)=-(AK(I,J,K)+AK(I-1,J,K))*
     1                   (X1(I,J,K)-X1(I-1,J,K))*R2DX2
            ENDDO
         ENDDO
      ENDDO


c$doacross local(k,j,i)
      DO K=2,KLES
         DO J=2,JMAX
            DO I=2,ILES
               Y3D(I,J,K)=-(AK(I,J,K)+AK(I,J-1,K))*
     1                      (X1(I,J,K)-X1(I,J-1,K))*R2DY2
            ENDDO
         ENDDO
      ENDDO


c$doacross local(k,j,i)
      DO K=2,KLES
         DO J=2,JLES
            DO I=2,ILES
               Y(I,J,K)=Y(I,J,K)+(Y1D(I,J,K)-Y1D(I+1,J,K))
     1                 +(Y3D(I,J,K)-Y3D(I,J+1,K))
            ENDDO
         ENDDO
      ENDDO

cccccccccccccccccccccc

c$doacross local(k,a4hk,j,i,ip,im,ipp)
      DO K=2,KLES
         a4hk=awks*a4h(k)
         do j=2,jles
            DO I=2,IL2
               IP=I+1
               IM=I-1
               IPP=I+2
               Y2D(I,J,K)=a4hk*((X1(IPP,J,K)-X1(IP,J,K))-2.*
     1            (X1(IP,J,K)-X1(I,J,K))+(X1(I,J,K)-X1(IM,J,K)))
            enddo
         enddo
      enddo


c$doacross local(k,a4hk,j)
      DO K=2,KLES
         a4hk=awks*a4h(k)
         do j=2,jles
            Y2D(ILES,J,K)=a4hk*((X1(3,J,K)-X1(IMAX,J,K))-2.*
     1         (X1(IMAX,J,K)-X1(ILES,J,K))+(X1(ILES,J,K)-X1(IL2,J,K)))
         enddo
      enddo


c$doacross local(k,j,i,im)
      DO K=2,KLES
         do j=2,jles
            DO I=3,ILES
               IM=I-1
               Y(I,J,K)=Y(I,J,K)+Y2D(IM,J,K)-Y2D(I,J,K)
            ENDDO
         enddo
      enddo


c$doacross local(k,j)
      DO K=2,KLES
         do j=2,jles
            Y(2,J,K)=Y(2,J,K)+Y2D(ILES,J,K)-Y2D(2,J,K)
         enddo
      enddo



c$doacross local(k,a4hk,j,jpp,i)
      DO K=2,KLES
         A4HK=awks*A4H(K)
         DO J=2,JL2
            JPP=J+2
            do i=2,iles
               Y1D(I,J,K)=A4HK*((X1(I,JPP,K)-X1(I,J+1,K))-2.*
     1            (X1(I,J+1,K)-X1(I,J,K))+(X1(I,J,K)-X1(I,J-1,K)))
            enddo
         enddo
      enddo


c$doacross local(k,a4hk,i)
      DO K=2,KLES
         A4HK=awks*A4H(K)
         do i=2,iles
            Y1D(I,JLES,K)=A4HK*((X1(I,3,K)-X1(I,JMAX,K))-2.*
     1         (X1(I,JMAX,K)-X1(I,JLES,K))+(X1(I,JLES,K)-X1(I,JL2,K)))
         enddo
      enddo

c$doacross local(k,j,jm,i)
      DO K=2,KLES
         DO J=3,JLES
            JM=J-1
            do i=2,iles
               Y(I,J,K)=Y(I,J,K)+Y1D(I,JM,K)-Y1D(I,J,K)
            ENDDO
         ENDDO
      enddo

c$doacross local(k,i)
      DO K=2,KLES
         do i=2,iles
            Y(I,2,K)=Y(I,2,K)+Y1D(I,jles,K)-Y1D(I,2,K)
         ENDDO
      enddo



c$doacross local(k,amkk,rfak,wbk_wbkp,kp,km,j,i)
         DO K=2,KLES
            amkk=a0000*AM(K)*RD4Z
            rfak=a_irf*rfa(k)
c            wbkp=wb(k+1)
c            wbk=wb(k)
            wbk_wbkp=wb(k+1)+wb(k)
            kp=k+1
            km=k-1
            DO J=2,JLES
               DO I=2,ILES
                  Y(I,J,K)=Y(I,J,K)-amkk*y03d(i,j,1)*wbk_wbkp
     1                    *(X(I,J,KP)-X(I,J,KM))-RFAK*X1(I,J,K)
               ENDDO
            ENDDO
         ENDDO


c$doacross local(k,j,i)
      DO K=2,KLES
         DO J=2,JLES
            DO I=2,ILES
               XX0(I,J,K)=Y(I,J,K)
               xx000(i,j,k)=y(i,j,k)-xx000(i,j,k)
            ENDDO
         ENDDO
      ENDDO



C
       IF (ISMG.EQ.5) THEN
          CALL FADVECT (X,Y,VTP)
c$doacross local(k,j,i)
          DO K=2,KLES
             DO J=2,JLES
                DO I=2,ILES
                   xx00(i,j,k)=y(i,j,k)-xx0(i,j,k)
                ENDDO
             ENDDO
         ENDDO
       ENDIF



      if (id .eq. 1 .or. idq .eq. 1) then
    
         do k=2,kles

c$doacross local(j,i)
               do j=2,jles
                  do i=2,iles
                     ibz(i,j,1)=1
                     ibz(i,j,2)=ics5(i,j,2)
                     ibz(i,j,3)=ics5(i,j,3)
                     ibz(i,j,4)=ics5(i,j,4)
                  enddo
               enddo

cc$doacross local(mt,a22,a33,j,i)
               do mt=1,4
                  a22=0.
                  a33=0.
                  do j=2,jles
                     do i=2,iles
                        a33=a33+xx000(I,J,K)*ibz(i,j,mt)
                        if(ismg.eq.5)a22=a22+xx00(i,j,k)*ibz(i,j,mt)
                     enddo
                  enddo
                  if (ismg.eq.1) sntd(k,mt,4)=sntd(k,mt,4)+a33*a20
                  if (ismg.eq.2) snqd(k,mt,4)=snqd(k,mt,4)+a33*a20
                  if (ismg.gt.2) snhd(k,mt,4)=snhd(k,mt,4)+a33*a20
                  if (ismg.eq.5) snhv(k,mt,4)=snhv(k,mt,4)+a22*a20
               enddo
            enddo


       do kc=1,3
       DO K=2,kles
c$doacross local(j,i,y200)
        do j=2,jles
          do i=2,iles
             y200=.01*(ww1(i,j,k)+wb(k)*y03d(i,j,1))
             ibz(i,j,1)=1
             ibz(i,j,2)=ics5(i,j,2)
             ibz(i,j,3)=ics5(i,j,3)
             ibz(i,j,4)=ics5(i,j,4)
        if (y200.lt.rby(kc)) then
           ibz(i,j,1)=0
           ibz(i,j,2)=0
           ibz(i,j,3)=0
           ibz(i,j,4)=0
        endif
        enddo
        enddo
cc$doacross local(mt,a22,a33,j,i)
               do mt=1,4
                  a22=0.
                  a33=0.
                  do j=2,jles
                     do i=2,iles
                        a33=a33+xx000(I,J,K)*ibz(i,j,mt)
                        if(ismg.eq.5)a22=a22+xx00(i,j,k)*ibz(i,j,mt)
                     enddo
                  enddo
                  if (ismg.eq.1) sntd(k,mt,kc)=sntd(k,mt,kc)+a33*a20
                  if (ismg.eq.2) snqd(k,mt,kc)=snqd(k,mt,kc)+a33*a20
                  if (ismg.gt.2) snhd(k,mt,kc)=snhd(k,mt,kc)+a33*a20
                  if (ismg.eq.5) snhv(k,mt,kc)=snhv(k,mt,kc)+a22*a20
               enddo
            enddo
         enddo


         do kc=5,7
            DO K=2,kles
c$doacross local(j,i,y200)
       do j=2,jles
          do i=2,iles
             y200=.01*(ww1(i,j,k)+wb(k)*y03d(i,j,1))
             ibz(i,j,1)=1
             ibz(i,j,2)=ics5(i,j,2)
             ibz(i,j,3)=ics5(i,j,3)
             ibz(i,j,4)=ics5(i,j,4)
        if (y200.gt.rby(kc)) then
           ibz(i,j,1)=0
           ibz(i,j,2)=0
           ibz(i,j,3)=0
           ibz(i,j,4)=0
        endif
        enddo
        enddo
cc$doacross local(mt,a22,a33,j,i)
               do mt=1,4
                  a22=0.
                  a33=0.
                  do j=2,jles
                     do i=2,iles
                        a33=a33+xx000(I,J,K)*ibz(i,j,mt)
                        if(ismg.eq.5)a22=a22+xx00(i,j,k)*ibz(i,j,mt)
                     enddo
                  enddo
            if (ismg.eq.1) sntd(k,mt,kc)=sntd(k,mt,kc)+a33*a20
            if (ismg.eq.2) snqd(k,mt,kc)=snqd(k,mt,kc)+a33*a20
            if (ismg.gt.2) snhd(k,mt,kc)=snhd(k,mt,kc)+a33*a20
            if (ismg.eq.5) snhv(k,mt,kc)=snhv(k,mt,kc)+a22*a20
         enddo
         enddo
         enddo

         endif

C
        CALL ADVECTN (X,X1)
C


c$doacross local(k,j,i)
      DO K=2,KLES
         DO J=2,JLES
            DO I=2,ILES
               X(I,J,K)=X(I,J,K)+Y(I,J,K)*DT
            ENDDO
         ENDDO
      ENDDO
c
c
      endif  !!sh (end of if (ijkadv=0) else (ijkadv=1))
c!!sh  (stop here 11/20/98)
c     ***********************************************************

      if(ismg .gt. 2) then

      do k=kt1,kt2
         y1(k)=0.
         y2(k)=0.
         y3(k)=rho(k)*dz/am(k)
      enddo

c$doacross local(k,j,i)
      do k=kt1,kt2
         do j=2,jles
            do i=2,iles
               y1(k)=y1(k)+max(x(i,j,k),0.0)
               y2(k)=y2(k)+max(-x(i,j,k),0.0)
            enddo
         enddo
      enddo

      do k=kt1,kt2
         sqa(k)=sqa(k)+y2(k)*dt/d22t*rijl2
      enddo

      a1=0.
      a2=0.

      do k=kt1,kt2
         a1=a1+y1(k)*y3(k)
         a2=a2+y2(k)*y3(k)
      enddo

      a0=0.
      if(a1.ne.0.) a0=(a1-a2)/a1

c$doacross local(k,j,i)
      do k=kt1,kt2
         do j=1,jmax
            do i=1,imax
               x(i,j,k)=max(x(i,j,k),0.0)*a0
            enddo
         enddo
      enddo

      endif

      d2t=d22t

      return
      end

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ADVECTN (X,X1)
C     ****   COMPUTE ADVECTION FOR VARIOUS SCALAR VARIABLES
      PARAMETER (NX=130,NY=130,NZ=43)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
      DIMENSION X(NX,NY,1),X1(NX,NY,1)
CCC   LOCAL VARIABLES
      COMMON /tmp1/ U1(NX,NY,NZ),V1(NX,NY,NZ),W1(NX,NY,NZ)
      COMMON/SAVEC/ ACONST
c     dan new
      save
CC    **************************************************************
c$doacross local(k,j,i)
      DO K=1,KMAX
      DO J=1,JMAX
      DO I=1,IMAX
	U1(I,J,K)=.5*(3.*UU1(I,J,K)-UMD(I,J,K))
	V1(I,J,K)=.5*(3.*VV1(I,J,K)-VMD(I,J,K))      
	W1(I,J,K)=.5*(3.*WW1(I,J,K)-WMD(I,J,K))
      ENDDO
      ENDDO
      ENDDO
c
	ACONST=0.0
       IF (ISMG.LE.2) THEN
	CALL TMAXADV (X,SMALL)
	 ACONST=ABS(SMALL*1.25)
c
c$doacross local(k,j,i)
	 DO K=1,KMAX
	    DO J=1,JMAX
	       DO I=1,IMAX
		  X(I,J,K)=X(I,J,K)+ACONST
	       ENDDO
	    ENDDO
	 ENDDO
      ENDIF

c$doacross local(k,j,i)
      DO K=1,KMAX
	 DO J=1,JMAX
	    DO I=1,IMAX
	       X1(I,J,K)=X(I,J,K)
	    ENDDO
	 ENDDO
      ENDDO

      CALL FADV (X,U1,V1,W1)
CC    NEED TO COMPUTE BOUNDARY CONDITION
       call bndop (x,x1)
CC    *********************************************
       CALL FADVUW (X,X1)
       CALL FADV (X,U1,V1,W1)
	call bndop (x,x1)
C
      IF(ISMG.LE.2) THEN
C
c$doacross local(k,j,i)
       DO K=1,KMAX
       DO J=1,JMAX
       DO I=1,IMAX
	X(I,J,K)=X(I,J,K)-ACONST
      ENDDO
      ENDDO
      ENDDO

c$doacross local(k,j,i)
       DO K=1,KMAX
       DO J=1,JMAX
       DO I=1,IMAX
	X1(I,J,K)=X1(I,J,K)-ACONST
      ENDDO
      ENDDO
      ENDDO

      ENDIF
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FADV (X,U,V,W)
C     ****   COMPUTE ADVECTION OF DIFFERENT TRACERS
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NM7=7*NM)
       common /setbound/ icheck
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      common/timestat/ ndt_stat,idq
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BT/ DT,D2T,RIJL2,DTS(16)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),
     1  TA1(NZ),QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),
     2  UB(NZ),VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/BCOR1/ DXR(NX,NY,1),DXXT(NX,NY),DXXT1(NX,NY),DYR(NX,NY,1),
     1   DYYT(NX,NY),DYYT1(NX,NY),DZR(NZ),DZZT(NZ),DZZT1(NZ)
      COMMON/B7/ SQAAQ(NZ),SQTT(NZ),SQVV(NZ),SQAAT(NZ),SQAK(NZ)
cc      COMMON/BW/ TRAHX(NX,NY,NZ),TRAHY(NX,NY,NZ),TRAV(NX,NY,NZ)
      COMMON/BBA/ XY1(NX,NY),XY2(NX,NY),XY3(NX,NY),XY4(NX,NY),Y1(NM),
     1   Y2(NM),Y3(NM),Y4(NM),Y5(NM),Y6(NM),Y7(NM),Y8(NM),Y9(NM),
     1   Y10(NM7)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/BSTS2/ snth(nz,4,7),sntv(nz,4,7),sntd(nz,4,7),SNTL(NZ,4,7),
     1 SNQH(NZ,4,7),snqv(nz,4,7),snqd(nz,4,7),SNQL(NZ,4,7),
     2 SNTL1(NZ,4,7),snhh(nz,4,7),snhv(nz,4,7),snhd(nz,4,7),
     3 SNQL1(NZ,4,7)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NX,NY),ICS5(NX,NY,4),
     1  IBZ(NX,NY,4)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
      common/bch1/ rby(7)
      common/TMP/ y1d(nx,ny,nz),y2d(nx,ny,nz),y3d(nx,ny,nz),
     1            y4d(nx,ny,nz),y03d(nx,ny,1)
      COMMON/BSAT1/ AAA(NX,NY,NZ)

C_TAO

      DIMENSION U(NX,NY,1),V(NX,NY,1),W(NX,NY,1),X(NX,NY,1)
C     LOCAL VARIABLES
C     ******   HORIZONTAL AND VERTICAL ADVECTION TERMS   ***************
c     dan new
      save

      a20=ndt_stat*rijl2

	 kbot=2
	 ktop=kles

       if(icheck.eq.1)then

      sum=0.
      do k=nz,1,-1
c$doacross local(j,i),reduction(sum)
	 do j=1,ny
	    do i=1,nx
	        sum=sum+abs(x(i,j,k))
	    enddo
	 enddo
	 if(sum.ne.0.)then
	     ktop=k+1
	     goto 2222
	 endif
      enddo
      ktop=0
2222   continue
       ktop=min(kles,ktop)


      sum=0.
      do k=1,ktop
c$doacross local(j,i),reduction(sum)
	 do j=1,ny
	    do i=1,nx
	        sum=sum+abs(x(i,j,k))
	    enddo
	 enddo
	 if(sum.ne.0.)then
	     kbot=k-1
	     goto 3333
	 endif
      enddo
      kbot=kles
3333   continue
       kbot=max(2,kbot)

       endif


c$doacross local(k,km,rho1k,j,jm,i,im,tmpu,tmpv,tmpw)
      do k=kbot,ktop+1
         km=k-1
         rho1k=rho1(k)
         do j=2,ny
            jm=j-1
            do i=2,nx
               im=i-1
               tmpu=x(i,j,k)
               tmpv=x(i,j,k)
               tmpw=x(i,j,k)
               if(u(i,j,k).gt.0.)tmpu=x(im,j,k)
               if(v(i,j,k).gt.0.)tmpv=x(i,jm,k)
               if(w(i,j,k).gt.0.)tmpw=x(i,j,km)
               y2d(i,j,k)=rho1k*w(i,j,k)*tmpw
               y3d(i,j,k)=u(i,j,k)*tmpu
               y4d(i,j,k)=v(i,j,k)*tmpv
            enddo
         enddo
      enddo

c$doacross local(k,kp,dzrk,j,jp,i)
      DO K=kbot,ktop
         KP=K+1
         dzrk=dzr(k)
	 DO J=2,JLES
            JP=J+1
	    DO I=2,ILES
              Y1D(I,J,K)=DXR(I,J,1)*(Y3D(I,J,K)-Y3D(I+1,J,K))+
     1                    DYR(I,J,1)*(Y4D(I,J,K)-Y4D(I,JP,K))
                AAA(I,J,K)=DZRK*(Y2D(I,J,K)-Y2D(I,J,KP))
	        X(I,J,K)=X(I,J,K)+(Y1D(I,J,K)+AAA(I,J,K))*DT 
            ENDDO
         ENDDO
      ENDDO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(id.eq.1 .or. idq .eq. 1) then


	    DO K=kbot,ktop

c$doacross local(j,i)
	       do j=2,jles
		  do i=2,iles
		     ibz(i,j,1)=1
		     ibz(i,j,2)=ics5(i,j,2)
		     ibz(i,j,3)=ics5(i,j,3)
		     ibz(i,j,4)=ics5(i,j,4)
		  enddo
	       enddo

c$doacross local(mt,a11,a22,j,i)
	       do mt=1,4
		  a11=0.
		  a22=0.
cc$doacross local(j,i),reduction(a11,a22)
		  do j=2,jles
		     do i=2,iles
		      a11=a11+y1d(i,j,k)*ibz(i,j,mt)
		      a22=a22+aaa(i,j,k)*ibz(i,j,mt)
		     enddo
		  enddo

		  if (ismg.eq.1) then
		     snth(k,mt,4)=snth(k,mt,4)+a11*a20
		     sntv(k,mt,4)=sntv(k,mt,4)+a22*a20
		  elseif (ismg.eq.2) then
		     snqh(k,mt,4)=snqh(k,mt,4)+a11*a20
		     snqv(k,mt,4)=snqv(k,mt,4)+a22*a20
		  elseif (ismg.gt.2) then
		     snhh(k,mt,4)=snhh(k,mt,4)+a11*a20
		     snhv(k,mt,4)=snhv(k,mt,4)+a22*a20
		  endif
	       enddo
	    enddo


	 do kc=1,3
	    DO K=kbot,ktop
c$doacross local(j,i,y100)
       do j=2,jles
	  do i=2,iles
	     y100=.01*(ww1(i,j,k)+wb(k)*y0(i,j))
	     ibz(i,j,1)=1
	     ibz(i,j,2)=ics5(i,j,2)
	     ibz(i,j,3)=ics5(i,j,3)
	     ibz(i,j,4)=ics5(i,j,4)
	if (y100.lt.rby(kc)) then
	   ibz(i,j,1)=0
	   ibz(i,j,2)=0
	   ibz(i,j,3)=0
			   ibz(i,j,4)=0
			endif
		  enddo
	       enddo

c$doacross local(mt,a11,a22,j,i)
	       do mt=1,4
		  a11=0.
		  a22=0.
cc$doacross local(j,i),reduction(a11,a22)
		  do j=2,jles
		     do i=2,iles
			if(ibz(i,j,mt).eq.1) then
			   a11=a11+y1d(i,j,k)
			   a22=a22+aaa(i,j,k)
			endif
		     enddo
		  enddo
		  if (ismg.eq.1) then
		     snth(k,mt,kc)=snth(k,mt,kc)+a11*a20
		     sntv(k,mt,kc)=sntv(k,mt,kc)+a22*a20
		  elseif (ismg.eq.2) then
		     snqh(k,mt,kc)=snqh(k,mt,kc)+a11*a20
		     snqv(k,mt,kc)=snqv(k,mt,kc)+a22*a20
		  elseif (ismg.gt.2) then
		     snhh(k,mt,kc)=snhh(k,mt,kc)+a11*a20
		     snhv(k,mt,kc)=snhv(k,mt,kc)+a22*a20
		  endif
	       enddo
	    enddo
	 ENDDO

	 do kc=5,7
	    DO K=kbot,ktop
c$doacross local(j,i,y100)
	       do j=2,jles
		  do i=2,iles
		     y100=.01*(ww1(i,j,k)+wb(k)*y0(i,j))
		     ibz(i,j,1)=1
		     ibz(i,j,2)=ics5(i,j,2)
		     ibz(i,j,3)=ics5(i,j,3)
		     ibz(i,j,4)=ics5(i,j,4)
			if (y100.gt.rby(kc)) then
			   ibz(i,j,1)=0
			   ibz(i,j,2)=0
			   ibz(i,j,3)=0
			   ibz(i,j,4)=0
			endif
		  enddo
	       enddo

c$doacross local(mt,a11,a22,j,i)
	       do mt=1,4
		  a11=0.
		  a22=0.
cc$doacross local(j,i),reduction(a11,a22)
		  do j=2,jles
		     do i=2,iles
			if(ibz(i,j,mt).eq.1) then
			   a11=a11+y1d(i,j,k)
			   a22=a22+aaa(i,j,k)
			endif
		     enddo
		  enddo
		  if (ismg.eq.1) then
		     snth(k,mt,kc)=snth(k,mt,kc)+a11*a20
		     sntv(k,mt,kc)=sntv(k,mt,kc)+a22*a20
		  elseif (ismg.eq.2) then
		     snqh(k,mt,kc)=snqh(k,mt,kc)+a11*a20
		     snqv(k,mt,kc)=snqv(k,mt,kc)+a22*a20
		  elseif (ismg.gt.2) then
		     snhh(k,mt,kc)=snhh(k,mt,kc)+a11*a20
		     snhv(k,mt,kc)=snhv(k,mt,kc)+a22*a20
		  endif
	       enddo
	    enddo
	 ENDDO

      endif

      RETURN
      END


      SUBROUTINE FADVUW (X,X1)
CC    ****   COMPUTE ANTI-DIFFUSIVE VELOCITIES   *******************
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/OPTION1/ ICOMP,NSML
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS(16)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/BBA/ XY1(NX,NY),XY2(NX,NY),XY3(NX,NY),XY4(NX,NY),Y1(NM),
     1   Y2(NM),Y3(NM),Y4(NM),Y5(NM),Y6(NM),Y7(NM),Y8(NM),Y9(NM),
     1   Y10(NM),Y11(NM),Y12(NM),Y13(NM),Y14(NM),Y15(NM),Y16(NM)
      COMMON/BCOR1/ DXR(NX,NY),DXXT(NX,NY,1),DXXT1(NX,NY),DYR(NX,NY),
     1   DYYT(NX,NY,1),DYYT1(NX,NY),DZR(NZ),DZZT(NZ),DZZT1(NZ)
      common/TMP/ UMOD(nx,ny,nz),VMOD(nx,ny,nz),y3d(nx,ny,nz),
     1            y4d(nx,ny,nz),y03d(nx,ny,1)
      DIMENSION X(NX,NY,1),X1(NX,NY,1)
      COMMON/TMP1/ U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ)
      COMMON/BSAT1/ WMOD(NX,NY,NZ)
CC    LOCAL VARIABLES
      real tmp1(nx,1,nz),tmp2(nx,1,nz)
      real umod2(nx,ny,nz),vmod2(nx,ny,nz),wmod2(nx,ny,nz)
      save

       EPS=1.E-10
       A18X=1./8.*DT*RDX
       A18Y=1./8.*DT*RDY
       A18Z=1./8.*DT*RDZ

c$doacross local(k,j,i)
       DO K=1,KMAX
	  DO J=1,JMAX
	     DO I=1,IMAX
		UMOD(I,J,K)=0.
		VMOD(I,J,K)=0.
		WMOD(I,J,K)=0.
	     enddo
	  enddo
       enddo


C   ***   U-COMPONENT WIND   *******************************************
c$doacross local(k,kp,km,kmm,ar1,ar1k,arp,arm,a00,a0x,a0y,j,jp,jm,jmm,i,
c$&              ip,im,imm,uiip,vjjp,y55,y44,y66,y144,y166,ximmjk)
      DO K=2,KLES
	 KP=K+1
         KM=K-1
	 KMM=K-2
	 AR1=RHO1(K)
	 AR1K=RHO(K)
	 ARP=RHO1(KP)
	 ARM=RHO(KM)
	 A00=A18Z*AM(K)*RRHO(K)
	 A0X=A18X*RRHO1(K)
	 A0Y=A18Y*RRHO1(K)

	 DO J=2,JLES
	    jp=j+1
	    jm=j-1
	    jmm=j-2
	    if(j.eq.2)jmm=jl2

	    DO I=2,ILES
	       IP=I+1
	       IM=I-1
	       imm=i-2
	       if(i.eq.2) imm=il2
               UIIP=U(I,J,K)+U(IP,J,K)
               VJJP=V(I,J,K)+V(I,JP,K)

	       Y55=X(I,J,K)+X(IM,J,K)

	       IF (Y55 .GE. EPS) THEN
		  Y44=X(I,J,KP)+X(IM,J,KP)
		  Y66=X(I,J,KM)+X(IM,J,KM)
		  Y144=X(I,JP,K)+X(IM,JP,K)
		  Y166=X(I,JM,K)+X(IM,JM,K)
                  XIMMJK=X(IMM,J,K)+X(IP,J,K) 

	       UMOD(I,J,K)=(ABS(U(I,J,K))-U(I,J,K)*U(I,J,K)*DXXT(I,J,1))
     1           *(X(I,J,K)-X(IM,J,K))/(Y55+EPS)
c     2           -A00*U(I,J,K)*(Y44-Y66)/(Y44+Y66+EPS)
c     3           *(AR1*(W(I,J,K)+W(IM,J,K))+ARP*(W(I,J,KP)+W(IM,J,KP)))
c     4           -A18Y*U(I,J,K)*(Y144-Y166)/(Y144+Y166+EPS)
c     5           *(VJJP+V(IM,J,K)+V(IM,JP,K))
c     6           -U(I,J,K)*(XIMMJK-X(IM,J,K)-X(I,J,K))
c     7           /(3.*(XIMMJK+X(IM,J,K)+X(I,J,K)+EPS))
	       ENDIF


	       Y55=X(I,J,K)+X(I,JM,K)

	       IF (Y55 .GE. EPS) THEN
		  Y44=X(I,J,KP)+X(I,JM,KP)
		  Y66=X(I,J,KM)+X(I,JM,KM)
		  Y144=X(IP,J,K)+X(IP,JM,K)
		  Y166=X(IM,J,K)+X(IM,JM,K)

	       VMOD(I,J,K)=(ABS(V(I,J,K))-V(I,J,K)*V(I,J,K)*DYYT(I,J,1))
     1            *(X(I,J,K)-X(I,JM,K))/(Y55+EPS)
c     2            -A00*V(I,J,K)*(Y44-Y66)/(Y44+Y66+EPS)
c     3            *(AR1*(W(I,J,K)+W(I,JM,K))+ARP*(W(I,J,KP)+W(I,JM,KP)))
c     4            -A18X*V(I,J,K)*(Y144-Y166)/(Y144+Y166+EPS)
c     5            *(UIIP+U(I,JM,K)+U(IP,JM,K))
c     6            -V(I,J,K)*(X(I,JMM,K)-X(I,JM,K)-X(I,J,K)+X(I,JP,K))
c     7            /(3.*(X(I,JMM,K)+X(I,JM,K)+X(I,J,K)+X(I,JP,K)+EPS))
	       ENDIF


 	       if(k.gt.2)then
	       Y55=X(I,J,K)+X(I,J,KM)
C
	       IF (Y55 .GE. EPS) THEN
		  Y44=X(IP,J,K)+X(IP,J,KM)
		  Y66=X(IM,J,K)+X(IM,J,KM)
		  Y144=X(I,JP,K)+X(I,JP,KM)
		  Y166=X(I,JM,K)+X(I,JM,KM)

		 WMOD(I,J,K)=(ABS(W(I,J,K))-W(I,J,K)*W(I,J,K)*DZZT(K))
     1           *(X(I,J,K)-X(I,J,KM))/(Y55+EPS)
c     2           -A0X*W(I,J,K)*(Y44-Y66)/(Y44+Y66+EPS)
c     3           *(AR1K*UIIP+ARM*(U(I,J,KM)+U(IP,J,KM)))
c     4           -A0Y*W(I,J,K)*(Y144-Y166)/(Y144+Y166+EPS)
c     5           *(AR1K*VJJP+ARM*(V(I,J,KM)+V(I,JP,KM)))
c     6           -W(I,J,K)*(X(I,J,KMM)-X(I,J,KM)-X(I,J,K)+X(I,J,KP))
c     7           /(3.*(X(I,J,KMM)+X(I,J,KM)+X(I,J,K)+X(I,J,KP)+EPS))
	       ENDIF
	      endif
	    ENDDO
	 ENDDO
      ENDDO

      call bndop (umod,vmod)
      call bndop (wmod,wmod)

c$doacross local(j,i)
      do j=1,jmax
	 do i=1,imax
	    wmod(i,j,kmax)=0.
	 enddo
      enddo

CCCCC   NON-OSCILLATORY OPTION (SMOLARKIEWICZ AND GRABOWSKI, 1990)

c$doacross local(k,j,i)
      DO K=1,KMAX
	 DO J=1,JMAX      
	    DO I=1,IMAX
	       U(I,J,K)=0.
	       V(I,J,K)=0.
	       W(I,J,K)=0.
	       UMOD2(I,J,K)=MIN(UMOD(I,J,K),0.)
	       VMOD2(I,J,K)=MIN(VMOD(I,J,K),0.)
	       WMOD2(I,J,K)=MIN(WMOD(I,J,K),0.)
   	       UMOD(I,J,K)=MAX(UMOD(I,J,K),0.)
	       VMOD(I,J,K)=MAX(VMOD(I,J,K),0.)
	       WMOD(I,J,K)=MAX(WMOD(I,J,K),0.)
	    ENDDO
	 ENDDO
      ENDDO


c$doacross local(k,kp,km,dttk,j,jp,jm,i,im,ip,y33,y44,y55,y111,y222,
c$&        dtti,dttj,xx1max,xx1min,y3m,y4m,u1,u2)
      DO K=2,KLES
	 KP=K+1
	 KM=K-1
	 DTTK=DZZT1(K)

	 DO J=2,JLES
	    jp=j+1
	    jm=j-1
	    Y3M=0.
	    Y4M=0.
	    Y55=X(iles,J,K)+X(IL2,J,K)
	    IF (Y55 .GE. EPS) THEN
	       Y111=MAX(X(IL2,J,K),X(iles,J,K),X(NX,J,K),X1(IL2,J,K),
     1                   X1(iles,J,K),X1(NX,J,K))
	       Y222=MIN(X(IL2,J,K),X(iles,J,K),X(NX,J,K),X1(IL2,J,K),
     1                   X1(iles,J,K),X1(NX,J,K))
	       DTTI=DXXT1(iles,J)
	       DTTJ=DYYT1(iles,J)
	 if(y111.ne.x(iles,j,k))then
	Y3M=(Y111-X(iles,J,K))/(DTTI*(UMOD(ILES,J,K)*X(IL2,J,K)
     1   -UMOD2(NX,J,K)*X(NX,J,K))+DTTJ*(VMOD(iles,J,K)*X(ILES,JM,K)
     3   -VMOD2(iles,JP,K)*X(ILES,JP,K))+DTTK*(WMOD(iles,J,K)
     4    *X(ILES,J,KM)-WMOD2(iles,J,KP)*X(ILES,J,KP))+EPS)
	    endif
	   if(x(iles,j,k).ne.y222)then
	Y4M=(X(iles,J,K)-Y222)/(DTTI*(UMOD(NX,J,K)*X(ILES,J,K)
     1  -UMOD2(iles,J,K)*X(ILES,J,K))+DTTJ*(VMOD(iles,JP,K)*X(ILES,J,K)
     3  -VMOD2(iles,J,K)*X(ILES,J,K))+DTTK*(WMOD(iles,J,KP)*X(ILES,J,K)
     5                 -WMOD2(iles,J,K)*X(ILES,J,K))+EPS)
	    endif
	    ENDIF

	    DO I=2,ILES
	       IM=I-1
	       IP=I+1
	       DTTI=DXXT1(I,J)
	       DTTJ=DYYT1(I,J)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
              if (j.eq.2)then
               tmp1(I,1,K)=0.
               tmp2(I,1,K)=0.
	       Y55=X(I,JLES,K)+X(I,JL2,K)
	       IF (Y55 .GE. EPS) THEN
		  Y111=MAX(X(I,JL2,K),X(I,JLES,K),X(I,NY,K),X1(I,JL2,K),
     1                     X1(I,JLES,K),X1(I,NY,K))
		  Y222=MIN(X(I,JL2,K),X(I,JLES,K),X(I,NY,K),X1(I,JL2,K),
     1                     X1(I,JLES,K),X1(I,NY,K))
		  DTTI=DXXT1(I,JLES)
		  DTTJ=DYYT1(I,JLES)
	        if(y111.ne.x(i,jles,k))then
	        tmp1(I,1,K)=(Y111-X(I,JLES,K))/(DTTJ*(VMOD(I,JLES,K)
     1              *X(I,JL2,K)-VMOD2(I,NY,K)*X(I,NY,K))
     2                       +DTTI*(UMOD(I,JLES,K)*X(IM,JLES,K)
     3                       -UMOD2(IP,JLES,K)*X(IP,JLES,K))
     4                       +DTTK*(WMOD(I,JLES,K)*X(I,JLES,KM)
     5                       -WMOD2(I,JLES,KP)*X(I,JLES,KP))+EPS)
	 endif
	   if(x(i,jles,k).ne.y222)then
	tmp2(I,1,K)=(X(I,JLES,K)-Y222)/(DTTJ*(VMOD(I,NY,K)
     1               *X(I,JLES,K)-VMOD2(I,JLES,K)*X(I,JLES,K))
     2                       +DTTI*(UMOD(IP,JLES,K)*X(I,JLES,K)
     3                       -UMOD2(I,JLES,K)*X(I,JLES,K))
     4                       +DTTK*(WMOD(I,JLES,KP)*X(I,JLES,K)
     5                       -WMOD2(I,JLES,K)*X(I,JLES,K))+EPS)
	   endif
           ENDIF
               endif

               XX1MAX=MAX(X(I,J,K),X1(I,J,K))
               XX1MIN=MIN(X(I,J,K),X1(I,J,K))
               U1=DTTI*(UMOD(I,J,K)*X(IM,J,K)-UMOD2(IP,J,K)*X(IP,J,K))
     1          +DTTJ*(VMOD(I,J,K)*X(I,JM,K)-VMOD2(I,JP,K)*X(I,JP,K))
     2          +DTTK*(WMOD(I,J,K)*X(I,J,KM)-WMOD2(I,J,KP)*X(I,J,KP))
     3          +EPS
               U2=DTTI*(UMOD(IP,J,K)*X(I,J,K)-UMOD2(I,J,K)*X(I,J,K))
     1          +DTTJ*(VMOD(I,JP,K)*X(I,J,K)-VMOD2(I,J,K)*X(I,J,K))
     2          +DTTK*(WMOD(I,J,KP)*X(I,J,K)-WMOD2(I,J,K)*X(I,J,K))
     3          +EPS
c ------------------------------
	       Y33=0.
	       Y44=0.
	       Y55=X(I,J,K)+X(IM,J,K)
	       IF (Y55 .GE. EPS) THEN
              Y111=MAX(X(IM,J,K),XX1MAX,X(IP,J,K),X1(IM,J,K),X1(IP,J,K))
	      Y222=MIN(X(IM,J,K),XX1MIN,X(IP,J,K),X1(IM,J,K),X1(IP,J,K))
                 Y33=(Y111-X(I,J,K))/U1
                 Y44=(X(I,J,K)-Y222)/U2
	          U(I,J,K)=MIN(1.,Y4M,Y33)*UMOD(I,J,K)+MIN(1.,Y44,Y3M) 
     1                    *UMOD2(I,J,K)
	       ENDIF
	       Y3M=Y33 
	       Y4M=Y44
c ------------------------------
               Y33=0.
               Y44=0.
	       Y55=X(I,J,K)+X(I,JM,K)
	       IF (Y55 .GE. EPS) THEN
              Y111=MAX(X(I,JM,K),XX1MAX,X(I,JP,K),X1(I,JM,K),X1(I,JP,K))
              Y222=MIN(X(I,JM,K),XX1MIN,X(I,JP,K),X1(I,JM,K),X1(I,JP,K))
	       Y33=(Y111-X(I,J,K))/U1
	       Y44=(X(I,J,K)-Y222)/U2
	          V(I,J,K)=MIN(1.,tmp2(I,1,K),Y33)*VMOD(I,J,K)
     1                    +MIN(1.,Y44,tmp1(I,1,K))*VMOD2(I,J,K)
               ENDIF
               tmp1(I,1,K)=Y33
               tmp2(I,1,K)=Y44
c ------------------------------ 
	       y3d(i,j,k)=0.
	       y4d(i,j,k)=0.
	       Y55=X(I,J,K)+X(I,J,KM)
	       IF (Y55 .GE. EPS) THEN
              Y111=MAX(X(I,J,KM),XX1MAX,X(I,J,KP),X1(I,J,KM),X1(I,J,KP))
              Y222=MIN(X(I,J,KM),XX1MIN,X(I,J,KP),X1(I,J,KM),X1(I,J,KP))
              Y3D(I,J,K)=(Y111-X(I,J,K))/U1
	      Y4D(I,J,K)=(X(I,J,K)-Y222)/U2
	       ENDIF
            ENDDO
         ENDDO
      ENDDO

c$doacross local(k,km,j,i,y55)
      DO K=3,KLES
         KM=K-1
         DO J=2,JLES
            DO I=2,ILES
	       Y55=X(I,J,K)+X(I,J,KM)
	       IF (Y55 .GE. EPS) THEN
	          W(I,J,K)=MIN(1.,Y4D(I,J,KM),Y3D(I,J,K))*WMOD(I,J,K)
     1                    +MIN(1.,Y4D(I,J,K),Y3D(I,J,KM))*WMOD2(I,J,K)
	       ENDIF
	    ENDDO
	 ENDDO
      ENDDO

      call bndop (u,v)
      call bndop (w,w)

c$doacross local(j,i)
      do j=1,jmax
	 do i=1,imax
	    w(i,j,kmax)=0.
	 enddo
      enddo

      RETURN
      END


      SUBROUTINE FADVECT (X,Y,W)
C     ****   COMPUTE ADVECTION OF TERMINAL VELOCITY
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NZ2=2*NZ,NZ18=18*NZ,NB=6*NX*NY) 
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/B5/ TB(NZ2),RHO1(1,1,NZ),RHO(NZ18),WBX(NX)
      COMMON/BCOR1/ DXR(NB),DZR(NZ),DZZT(NZ),DZZT1(NZ)
      DIMENSION X(NX,NY,1),Y(NX,NY,1),W(NX,NY,1)
      save
C     ******   VERTICAL ADVECTION FOR TERMINAL VELOCITY   *************

c$doacross local(k,kp,dzrk,j,i)
      DO K=2,KLES
	 KP=K+1
	 DZRK=DZR(K)
	 DO J=2,JLES
	    DO I=2,ILES
	       Y(I,J,K)=Y(I,J,K)+DZRK*(RHO1(1,1,KP)*W(I,J,KP)*X(I,J,KP)
     1                 -RHO1(1,1,K)*W(I,J,K)*X(I,J,K))
	    ENDDO
	 ENDDO   
      ENDDO
      RETURN
      END
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SATICERH (FV)
C     (R&H)  COMPUTE ICE PHASE MICROPHYSICS AND SATURATION PROCESSES
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880)
      PARAMETER (NZ2=2*NZ,NZ3=3*NZ,NZ4=4*NZ,NM2=2*NM,NM11=11*NM)
      PARAMETER (NXY=NX*NY,NB=NX*NY*(NZ-29),NB2=NX*NY*(NZ-25),NT2=2*NT)
      PARAMETER (NB3=NZ+NX,NB4=5*NT,NB5=NX*NY*(NZ-29))
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      common/timestat/ ndt_stat,idq
C

cccshie by tao 5/3/01, shie 11/16/01 3d
      common/iice/ new_ice_sat

      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1  NRAN,KT1,KT2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/B3CS/ AG,BG,AS,BS,AW,BW,BGH,BGQ,BSH,BSQ,BWH,BWQ
      COMMON/SIZE/ TNW,TNS,TNG,ROQS,ROQG,ROQR
      COMMON/RTERV/ ZRC,ZGC,ZSC,VR0,VR1,VR2,VR3,VGC,VSC
      COMMON/RSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,RN2,BND2,RN3,RN4,
     1  RN5,RN50,RN51,RN52,RN53,RN6,RN60,RN61,RN62,RN63,RN7,RN8,RN9,
     2  RN10,RN101,RN102,RN10A,RN10B,RN10C,RN11,RN12,RN12A(31),
     3  RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN171,RN172,RN17A,RN17B,
     4  RN17C,RN18,RN18A,RN19,RN191,RN192,RN19A,RN20,RN20A,RN20B,RN30,
     5  RN30A,RN21,BND21,RN22,RN23,RN231,RN232,RN25,RN25A(31),RN31,BETA,
     6  RN32,RN33,RN331,RN332,RN34,RN35
C
      COMMON/B1TQ/ DPT(NXY,NZ),DQV(NXY,NZ)
      COMMON/B1CR/ QCL(NXY,NZ),QRN(NXY,NZ)
      COMMON/B1IG/ QCI(NXY,NZ),QCG(NXY,NZ)
      COMMON/B2TQ/ DPT1(NXY,NZ),DQV1(NXY,NZ)
      COMMON/B2CR/ QCL1(NXY,NZ),QRN1(NXY,NZ)
      COMMON/B2IG/ QCI1(NXY,NZ),QCG1(NXY,NZ)
      COMMON/B1S/ QCS(NXY,NZ)
      COMMON/B2S/ QCS1(NXY,NZ)
      COMMON/B4WP/ WW1(NXY,NZ)
      COMMON/SLWAVE/ RSW(NXY,NZ),RLW(NXY,NZ)
C
c      COMMON/BW/ TRAHX(NX,NY,NZ),TRAHY(NX,NY,NZ),TRAV(NX,NY,NZ)
C
      COMMON/BSAT1/ PT(NXY),QV(NXY),QC(NXY),QR(NXY),QI(NXY),QS(NXY),
     1 QG(NXY),TAIR(NXY),TAIRC(NXY),RTAIR(NXY),DEP(NXY),DD(NXY),
     2 DD1(NXY),QVS(NXY),DM(NXY),RQ(NXY),RSUB1(NXY),COL(NXY),CND(NXY),
     3 ERN(NXY),DLT1(NXY),DLT2(NXY),DLT3(NXY),DLT4(NXY),ZR(NXY),VR(NXY),
     4 ZS(NXY),VS(NXY),DBZ(NXY),DDA(NB)
      COMMON/BADV/ VG(NXY),ZG(NXY),PS(NXY),PG(NXY),PRN(NXY),PSN(NXY),
     1 PWACS(NXY),WGACR(NXY),PIDEP(NXY),PINT(NXY),QSI(NXY),SSI(NXY),
     2 ESI(NXY),ESW(NXY),QSW(NXY),PR(NXY),SSW(NXY),PIHOM(NXY),PIDW(NXY),
     3 PIMLT(NXY),PSAUT(NXY),QRACS(NXY),PSACI(NXY),PSACW(NXY),
     4 QSACW(NXY),PRACI(NXY),PMLTS(NXY),PMLTG(NXY),ASSS(NXY),DDE(NB5)
      COMMON/BSAT/ PRAUT(NXY),PRACW(NXY),PSFW(NXY),PSFI(NXY),DGACS(NXY),
     1 DGACW(NXY),DGACI(NXY),DGACR(NXY),PGACS(NXY),WGACS(NXY),
     2 QGACW(NXY),WGACI(NXY),QGACR(NXY),PGWET(NXY),PGAUT(NXY),
     3 PRACS(NXY),PSACR(NXY),QSACR(NXY),PGFR(NXY),PSMLT(NXY),PGMLT(NXY),
     4 PSDEP(NXY),PGDEP(NXY),PIACR(NXY),Y5(NXY),DDB(NB2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),ZX(NZ4),AM(NZ),ZQ(NZ3),WB(NZ),ZW(NZ2),RRHO(NZ),WBX(NB3)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      common/b66b/ s_dep(nz),s_sub(nz),s_qrs(nz),s_qrl(nz),s_mel(nz),
     1   s_frz(nz)
      COMMON/BRH1/ SRRO(NZ),QRRO(NZ),SQC(NZ),SQR(NZ),SQI(NZ),SQS(NZ),
     1   SQG(NZ),STQC(NZ),STQR(NZ),STQI(NZ),STQS(NZ),STQG(NZ),TTTD(NB4)
      COMMON/BSTS/ THOM(NZ,4,7),TDW(NZ,4,7),TMLT(NZ,4,7),SAUT(NZ,4,7),
     1 SACI(NZ,4,7),SACW(NZ,4,7),RACI(NZ,4,7),TACR(NZ,4,7),RAUT(NZ,4,7),
     2 RACW(NZ,4,7),SFW(NZ,4,7),SFI(NZ,4,7),GACS(NZ,4,7),GACW(NZ,4,7),
     3 GACI(NZ,4,7),GACR(NZ,4,7),GWET(NZ,4,7),GAUT(NZ,4,7),RACS(NZ,4,7),
     4 SACR(NZ,4,7),GFR(NZ,4,7),SMLT(NZ,4,7),GMLT(NZ,4,7),SDEP(NZ,4,7),
     5 SSUB(NZ,4,7),GSUB(NZ,4,7),PERN(NZ,4,7),D3RI(NZ,4,7),D3IR(NZ,4,7),
     6 D2SR(NZ,4,7),D2RS(NZ,4,7),GDRY(NZ,4,7),COC(NZ,4,7),COE(NZ,4,7),
     7 SMF0(NZ,4,7),QC0(NZ,4,7),QR0(NZ,4,7),QI0(NZ,4,7),QS0(NZ,4,7),
     8 QG0(NZ,4,7),SQC0(NZ,4,7),SQR0(NZ,4,7),SQI0(NZ,4,7),SQS0(NZ,4,7),
     9 SQG0(NZ,4,7),ERNS(NZ,4,7),WGRS(NZ,4,7),QSWS(NZ,4,7),TB00(NZ,4),
     1 QB00(NZ,4)
      COMMON/BSTS1/ TUT1(NZ,4,7),TUT2(NZ,4,7),TVT1(NZ,4,7),TVT2(NZ,4,7),
     1 TSTF(NZ,4,7),TSTF1(NZ,4,7),TSTF2(NZ,4,7),TSQF(NZ,4,7),
     2 QQQ(NZ,4,7),TSQF1(NZ,4,7),TSQF2(NZ,4,7),TSQQ(NZ,4,7),
     3 TSQQ1(NZ,4,7)
      COMMON/BSTS3/ QV0(NZ,4,7),TT0(NZ,4,7),SQV0(NZ,4,7),STT0(NZ,4,7),
     1 SGPT(NZ,4,7),SSGPT(NZ,4,7),SNQHD(NZ,4,7),SNQVD(NZ,4,7),
     2 Q1T(NZ,4,7),SNHDH(NZ,4,7),SQHDT(NZ,4,7),SQVDT(NZ,4,7)
      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7),SQTDT(NZ,4,7),SQHL(NZ,4,7)

      common/bsts40/ fcld(nz,4,7)

      COMMON/BCS/ S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     1 S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     2 S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     3 SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     4 SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     5 SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCSS/ SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     1 SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     2 SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     3 SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     4 SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     5 SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)

      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NXY),ICS5(NXY,4),
     1  IBZ(NXY,4)
      COMMON/BI/ IT(NXY),ICS(NXY,4)
      COMMON/RSTAT/ CSTTT(NXY),CSTT(NXY)
      COMMON/BLS/ Y0(NXY),TS0NEW(NX,NY),QSS0NEW(NX,NY)
      COMMON/BBA/ Y1(NXY),Y2(NXY),Y3(NXY),Y4(NXY),Y7(NM2),B0(NM),B1(NM),
     1   B2(NM),Y6(NM11)


c    Q Budget for each point

C_TAO
c     COMMON/Q_BUGT/ Q1_G_H(NX,NY,NZ),Q1_G_V(NX,NY,NZ),
c    1               Q1A_G_H(NX,NY,NZ),Q1A_G_V(NX,NY,NZ),
c    2               Q1_D_H(NX,NY,NZ),Q1_D_V(NX,NY,NZ),
c    3               Q1A_D_H(NX,NY,NZ),Q1A_D_V(NX,NY,NZ),
c    4               Q2_G_H(NX,NY,NZ),Q2_G_V(NX,NY,NZ),
c    5               Q2A_G_H(NX,NY,NZ),Q2A_G_V(NX,NY,NZ),
c    6               Q2_D_H(NX,NY,NZ),Q2_D_V(NX,NY,NZ),
c    7               Q2A_D_H(NX,NY,NZ),Q2A_D_V(NX,NY,NZ),
c    8               Q1_HYD(NXY,NZ),Q2_HYD(NXY,NZ),
c    9               Q1A_HYD(NXY,NZ),Q2A_HYD(NXY,NZ),
c    9               Q1_RAD(NX,NY,NZ),Q1A_RAD(NX,NY,NZ),
c    9               IBUDSEC,RBUD

C_TAO



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION FV(1),RBY(7),AA1(31),AA2(31)
c     real tttbud(NXY),qqqbud(NXY)
      DATA AA1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,
     1   .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,
     2   .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,
     3   .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,
     4   .6513E-6,.5956E-6,.5333E-6,.4834E-6/
      DATA AA2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
      DATA RBY/2.,1.,0.,0.,0.,-.5,-1./
      save
CC    ***   THREE CLASSES OF ICE-PHASE   *******************************


      D22T=D2T

      IF(IJKADV .EQ. 1) THEN
	 D2T=DT
      ELSE
	 D2T=D2T
      ENDIF
c
      IJLES=NX*(NY-1)-1
      ISTART=NX+2
c
      ijles=max(nx*(ny-1), ny*(nx-1))-1
      istart=min(nx, ny)+2
c
      CMIN=1.E-20
      CMIN1=1.E-12
      CMIN2=1.E-20
c      CMAX2=1.E20
      UCOR=3071.29/TNW**.75                                            
      UCOS=687.97*ROQS**.25/TNS**.75
      UCOG=687.97*ROQG**.25/TNG**.75
      UWET=4.464**.95

c$doacross local(ij)
      DO IJ=1,NXY
	 IT(IJ)=1
      ENDDO

      FT=DT/D2T
      RFT=RIJL2*FT
      A0=.5*ISTATMIN*rijl2
      RT0=1./(T0-T00)
      BS3=BS+3.
      BG3=BG+3.
      BSH5=2.5+BSH
      BGH5=2.5+BGH
      BS6=6.+BS
      BETAH=.5*BETA
      RDT=1./D2T
      R10T=RN10*D2T
      R11T=RN11*D2T
      R19T=RN19*D2T
      R19AT=RN19A*D2T
      R20T=RN20*D2T
      R23T=RN23*D2T
	R25A=RN25
      R30T=RN30*D2T
      R33T=RN33*D2T

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    RH-TYPE ICE SCHEME

      DO 1000 K=KT1,KT2
c         if (ijkadv .eq. 1) then
c           tb0=ta(k)
c           qb0=qa(k)
c         else
	 tb0=ta1(k)
	 qb0=qa1(k)
c         endif
c         P00=P0(K)
	 RP0=3.799052E3/P0(K)
	 PI0=PI(K)
	 PIR=1./(PI(K))
CRH       PR0=1./P0(K)
	 R00=RHO(K)
c         R0S=SQRT(RHO(K))
	 RR0=RRHO(K)
	 RRS=SRRO(K)
	 RRQ=QRRO(K)
	 FV0=FV(K)
	 FVS=SQRT(FV(K))
	 CP409=C409*PI0
	 CV409=C409*AVC
	 CP580=C580*PI0
	 CS580=C580*ASC
CRH       ALVR=R00*ALV
	 AFCP=AFC*PIR
	 AVCP=AVC*PIR
	 ASCP=ASC*PIR
	 ZRR=1.E5*ZRC*RRQ
	 ZSR=1.E5*ZSC*RRQ
	 ZGR=1.E5*ZGC*RRQ
	 VSCF=VSC*FV0
	 VGCF=VGC*FV0
CS        VGCF=VGC*RRS
	 R1R=RN1*RR0
	 R3F=RN3*FV0
	 R4F=RN4*FV0
	 R5F=RN5*FV0
	 R6F=RN6*FV0
	 R7RF=RN7*RR0*FV0
	 R8RF=RN8*RR0*FV0
	 R9RF=RN9*RR0*FV0
	 R101R=RN101*RR0
	 R102RF=RN102*RRS*FVS
	 R12R=RN12*R00
	 R14F=RN14*FV0
	 R15F=RN15*FV0
	 R16RF=RN16*RR0*FV0
	 R18R=RN18*RR0
	 R191R=RN191*RR0
	 R192RF=RN192*RRS*FVS
	 R22F=RN22*FV0
	 R231R=RN231*RR0
	 R232RF=RN232*RRS*FVS
	 R25RT=RN25*RR0*D2T
	 R31R=RN31*RR0
	 R32RT=RN32*D2T*RRS
	 R331R=RN331*RR0
	 R332RF=RN332*RRS*FVS
	 R34F=RN34*FV0
	 SCC=0.
	 SEE=0.

c$doacross local(ij)
	 DO 125 IJ=ISTART,IJLES

c            tttbud(ij)=dpt(ij,k)
c            qqqbud(ij)=dqv(ij,k)

	    PT(IJ)=DPT(IJ,K)
	    QV(IJ)=DQV(IJ,K)
	    QC(IJ)=QCL(IJ,K)
	    QR(IJ)=QRN(IJ,K)
	    QI(IJ)=QCI(IJ,K)
	    QS(IJ)=QCS(IJ,K)
	    QG(IJ)=QCG(IJ,K)
C            IF (QV(IJ)+QB0 .LE. 0.) QV(IJ)=-QB0
	    IF (QC(IJ) .LE. CMIN) QC(IJ)=0.0
	    IF (QR(IJ) .LE. CMIN) QR(IJ)=0.0
	    IF (QI(IJ) .LE. CMIN) QI(IJ)=0.0
	    IF (QS(IJ) .LE. CMIN) QS(IJ)=0.0
	    IF (QG(IJ) .LE. CMIN) QG(IJ)=0.0
C            XX0(IJ)=PT(IJ)
C            XX00(IJ)=QV(IJ)
	    TAIR(IJ)=(PT(IJ)+TB0)*PI0
	    TAIRC(IJ)=TAIR(IJ)-T0
	    ZR(IJ)=ZRR
	    VR(IJ)=0.0
	    ZS(IJ)=ZSR
	    VS(IJ)=0.0
	    ZG(IJ)=ZGR
	    VG(IJ)=0.0

	    IF (QR(IJ) .Gt. CMIN) THEN
	       DD(IJ)=R00*QR(IJ)
	       Y1(IJ)=SQRT(DD(IJ))
	       Y2(IJ)=SQRT(Y1(IJ))
	       ZR(IJ)=ZRC/Y2(IJ)
	       VR(IJ)=FV0*(VR0+VR1*Y2(IJ)+VR2*Y1(IJ)+VR3*Y1(IJ)*Y2(IJ))
	       VR(IJ)=MAX(VR(IJ), 0.0)
	    ENDIF

	    IF (QS(IJ) .Gt. CMIN) THEN
	       DD(IJ)=R00*QS(IJ)
	       Y1(IJ)=DD(IJ)**.25
	       ZS(IJ)=ZSC/Y1(IJ)
	       VS(IJ)=MAX(VSCF*DD(IJ)**BSQ, 0.)
	    ENDIF

	    IF (QG(IJ) .Gt. CMIN) THEN
	       DD(IJ)=R00*QG(IJ)
	       Y1(IJ)=DD(IJ)**.25
	       ZG(IJ)=ZGC/Y1(IJ)
	       VG(IJ)=max(VGCF*DD(IJ)**BGQ, 0.0)
	    ENDIF

	    IF (QR(IJ) .LE. CMIN1) VR(IJ)=0.0
	    IF (QS(IJ) .LE. CMIN1) VS(IJ)=0.0
	    IF (QG(IJ) .LE. CMIN1) VG(IJ)=0.0
  125    CONTINUE

C*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
C*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
C*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
C* 34 * PWACS : COLLECTION OF QS BY QC                            **34**
C*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
C*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**

c$doacross local(ij,rn1,bnd1)
      DO 150 IJ=ISTART,IJLES
       PSAUT(IJ)=0.0
       PSACI(IJ)=0.0
       PRACI(IJ)=0.0
       PIACR(IJ)=0.0
       PSACW(IJ)=0.0
       PWACS(IJ)=0.0
       QSACW(IJ)=0.0
       IF(TAIR(IJ).LT.T0) THEN
c	  Y1(IJ)=RDT*(QI(IJ)-R1R*EXP(BETA*TAIRC(IJ)))
c	 PSAUT(IJ)=MAX(Y1(IJ),0.0)
            rn1=1.e-3
            bnd1=6.e-4
            esi(ij)=exp(.025*tairc(ij))
           psaut(ij)=max(rn1*esi(ij)*(qi(ij)-bnd1) ,0.0)
       ENDIF

       IF(TAIR(IJ).LT.T0) THEN
	PSACI(IJ)=R3F*QI(IJ)/ZS(IJ)**BS3
	PSACW(IJ)=R4F*QC(IJ)/ZS(IJ)**BS3
	PWACS(IJ)=R34F*QC(IJ)/ZS(IJ)**BS6
	  Y1(IJ)=1./ZR(IJ)
	 Y2(IJ)=Y1(IJ)*Y1(IJ)
	 Y3(IJ)=Y1(IJ)*Y2(IJ)
	 DD(IJ)=R5F*QI(IJ)*Y3(IJ)*(RN50+RN51*Y1(IJ)+RN52*Y2(IJ)
     1                            +RN53*Y3(IJ))
	PRACI(IJ)=MAX(DD(IJ),0.0)
	   Y4(IJ)=Y3(IJ)*Y3(IJ)
	 DD1(IJ)=R6F*QI(IJ)*Y4(IJ)*(RN60+RN61*Y1(IJ)+RN62*Y2(IJ)
     1                             +RN63*Y3(IJ))
	PIACR(IJ)=MAX(DD1(IJ),0.0)
       ELSE
	QSACW(IJ)=R4F*QC(IJ)/ZS(IJ)**BS3
       ENDIF

C* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
C* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
	 PRAUT(IJ)=MAX(RN21*(QC(IJ)-BND21),0.0)
	  Y1(IJ)=1./ZR(IJ)
	 Y2(IJ)=Y1(IJ)*Y1(IJ)
	 Y3(IJ)=Y1(IJ)*Y2(IJ)
	 Y4(IJ)=R22F*QC(IJ)*Y3(IJ)*(RN50+RN51*Y1(IJ)+RN52*Y2(IJ)
     1                             +RN53*Y3(IJ))
	PRACW(IJ)=MAX(Y4(IJ),0.0)
C* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
C* 13 * PSFI : BERGERON PROCESSES FOR QS                          **13**

	  PSFW(IJ)=0.0
	  PSFI(IJ)=0.0
	   IF(TAIR(IJ).LT.T0) THEN
	     Y1(IJ)=MAX( MIN(TAIRC(IJ), -1.), -31.)
	     IT(IJ)=INT(ABS(Y1(IJ)))
	     Y1(IJ)=RN12A(IT(IJ))
	     Y2(IJ)=RN12B(IT(IJ))
	     Y3(IJ)=RN13(IT(IJ))
	    PSFW(IJ)=MAX(D2T*Y1(IJ)*(Y2(IJ)+R12R*QC(IJ))*QI(IJ),0.0)
	    PSFI(IJ)=Y3(IJ)*QI(IJ)
	   ENDIF
CTTT***** QG=QG+MIN(PGDRY,PGWET)
C*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
C* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
C* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
C* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
	 Y1(IJ)=ABS( VG(IJ)-VS(IJ) )
	 Y2(IJ)=ZS(IJ)*ZG(IJ)
	 Y3(IJ)=5./Y2(IJ)
	 Y4(IJ)=.08*Y3(IJ)*Y3(IJ)
	 Y5(IJ)=.05*Y3(IJ)*Y4(IJ)
	 Y2(IJ)=Y1(IJ)*(Y3(IJ)/ZS(IJ)**5+Y4(IJ)/ZS(IJ)**3+Y5(IJ)
     1                  /ZS(IJ))

	PGACS(IJ)=R9RF*Y2(IJ)
	DGACS(IJ)=PGACS(IJ)
CRH     WGACS(IJ)=10.*R9RF*Y2(IJ)
	WGACS(IJ)=0.0
	 Y1(IJ)=1./ZG(IJ)**BG3
	DGACW(IJ)=R14F*QC(IJ)*Y1(IJ)
	QGACW(IJ)=DGACW(IJ)
	DGACI(IJ)=R15F*QI(IJ)*Y1(IJ)
CRH     WGACI(IJ)=R15AF*QI(IJ)*Y1(IJ)
	WGACI(IJ)=0.0

	 Y1(IJ)=ABS( VG(IJ)-VR(IJ) )
	 Y2(IJ)=ZR(IJ)*ZG(IJ)
	 Y3(IJ)=5./Y2(IJ)
	 Y4(IJ)=.08*Y3(IJ)*Y3(IJ)
	 Y5(IJ)=.05*Y3(IJ)*Y4(IJ)
	 DD(IJ)=R16RF*Y1(IJ)*(Y3(IJ)/ZR(IJ)**5+Y4(IJ)/ZR(IJ)**3+Y5(IJ)
     1                       /ZR(IJ))
	DGACR(IJ)=MAX(DD(IJ),0.0)
	QGACR(IJ)=DGACR(IJ)

	 IF (TAIR(IJ) .GE. T0) THEN
	  DGACS(IJ)=0.0
CRH       WGACS(IJ)=0.0
	  DGACW(IJ)=0.0
	  DGACI(IJ)=0.0
CRH       WGACI(IJ)=0.0
	  DGACR(IJ)=0.0
	 ELSE
	  PGACS(IJ)=0.0
	  QGACW(IJ)=0.0
	  QGACR(IJ)=0.0
	 ENDIF
  150  CONTINUE

C*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
C* 17 * PGWET : WET GROWTH OF QG                                  **17**
CRH       PGWET(IJ)=0.0
CRH       IF (TAIR(IJ) .LT. T0) THEN
CRH         Y1(IJ)=1./(ALF+RN17C*TAIRC(IJ))
CRH         Y2(IJ)=RP0-(QV(IJ)+QB0)
CRH         Y3(IJ)=.78/ZG(IJ)**2+R17ARF/ZG(IJ)**BGH5
CRH         Y4(IJ)=RN171*Y2(IJ)-R172R*TAIRC(IJ)
CRH         DD(IJ)=Y1(IJ)*(Y4(IJ)*Y3(IJ)+(WGACI(IJ)
CRH  1                               +WGACS(IJ))*(ALF+RN17B*TAIRC(IJ)))
CRH        PGWET(IJ)=MAX(DD(IJ), 0.0)
CRH       ENDIF
C******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
CRH     WGACR(IJ)=PGWET(IJ)-DGACW(IJ)-WGACI(IJ)-WGACS(IJ)
CRH      Y2(IJ)=DGACW(IJ)+DGACI(IJ)+DGACR(IJ)+DGACS(IJ)
CRH      IF (PGWET(IJ) .GE. Y2(IJ)) THEN
CRH       WGACR(IJ)=0.0
CRH       WGACI(IJ)=0.0
CRH       WGACS(IJ)=0.0
CRH      ELSE
CRH       DGACR(IJ)=0.0
CRH       DGACI(IJ)=0.0
CRH       DGACS(IJ)=0.0
CRH      ENDIF
C*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
C* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
C* 17 * PGWET : WET GROWTH OF QG                                  **17**
C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
C********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************

c$doacross local(ij)
      DO 175 IJ=ISTART,IJLES
	PGWET(IJ)=0.0
	   Y1(IJ)=QC(IJ)/D2T
	  PSACW(IJ)=MIN(Y1(IJ), PSACW(IJ))
	  PRAUT(IJ)=MIN(Y1(IJ), PRAUT(IJ))
	  PRACW(IJ)=MIN(Y1(IJ), PRACW(IJ))
	  PSFW(IJ)= MIN(Y1(IJ), PSFW(IJ))
	  DGACW(IJ)=MIN(Y1(IJ), DGACW(IJ))
	  QSACW(IJ)=MIN(Y1(IJ), QSACW(IJ))
	  QGACW(IJ)=MIN(Y1(IJ), QGACW(IJ))

	Y1(IJ)=D2T*(PSACW(IJ)+PRAUT(IJ)+PRACW(IJ)+PSFW(IJ)+DGACW(IJ)
     1             +QSACW(IJ)+QGACW(IJ))

	QC(IJ)=QC(IJ)-Y1(IJ)
c
	if (qc(ij) .lt. 0.0) then
	   y2(ij)=1.
	    if (y1(ij) .ne. 0.) y2(ij)=qc(ij)/y1(ij)+1.
	   psacw(ij)=psacw(ij)*y2(ij)
	   praut(ij)=praut(ij)*y2(ij)
	   pracw(ij)=pracw(ij)*y2(ij)
	   psfw(ij)=psfw(ij)*y2(ij)
	   dgacw(ij)=dgacw(ij)*y2(ij)
	   qsacw(ij)=qsacw(ij)*y2(ij)
	   qgacw(ij)=qgacw(ij)*y2(ij)
	   qc(ij)=0.0
	 endif
c
	    Y1(IJ)=QI(IJ)/D2T
	   PSAUT(IJ)=MIN(Y1(IJ), PSAUT(IJ))
	   PSACI(IJ)=MIN(Y1(IJ), PSACI(IJ))
	   PRACI(IJ)=MIN(Y1(IJ), PRACI(IJ))
	   PSFI(IJ)= MIN(Y1(IJ), PSFI(IJ))
	   DGACI(IJ)=MIN(Y1(IJ), DGACI(IJ))
	   WGACI(IJ)=MIN(Y1(IJ), WGACI(IJ))

	Y1(IJ)=D2T*(PSAUT(IJ)+PSACI(IJ)+PRACI(IJ)+PSFI(IJ)+DGACI(IJ)
     1             +WGACI(IJ))

	QI(IJ)=QI(IJ)-Y1(IJ)
c
	 if (qi(ij) .lt. 0.0) then
	   y2(ij)=1.
	    if (y1(ij) .ne. 0.0) y2(ij)=qi(ij)/y1(ij)+1.
	   psaut(ij)=psaut(ij)*y2(ij)
	   psaci(ij)=psaci(ij)*y2(ij)
	   praci(ij)=praci(ij)*y2(ij)
	   psfi(ij)=psfi(ij)*y2(ij)
	   dgaci(ij)=dgaci(ij)*y2(ij)
	   wgaci(ij)=wgaci(ij)*y2(ij)
	   qi(ij)=0.0
	 endif
c
	 WGACR(IJ)=QGACR(IJ)+QGACW(IJ)
	DLT3(IJ)=0.0
	 IF (QR(IJ) .LT. 1.E-4) DLT3(IJ)=1.
	DLT4(IJ)=1.
	 IF (QC(IJ) .GT. 5.E-4) DLT4(IJ)=0.0
	 IF (QS(IJ) .LE. 1.E-4) DLT4(IJ)=1.
	  IF (TAIR(IJ) .GE. T0) THEN
	   DLT3(IJ)=0.0
	   DLT4(IJ)=0.0
	  ENDIF
	PR(IJ)=D2T*(QSACW(IJ)+PRAUT(IJ)+PRACW(IJ)+WGACR(IJ)-QGACR(IJ))
	PS(IJ)=D2T*(PSAUT(IJ)+PSACI(IJ)+DLT4(IJ)*PSACW(IJ)+PSFW(IJ)
     1             +PSFI(IJ)+DLT3(IJ)*PRACI(IJ))
	PG(IJ)=D2T*((1.-DLT3(IJ))*PRACI(IJ)+DGACI(IJ)+WGACI(IJ)
     1             +DGACW(IJ)+(1.-DLT4(IJ))*PSACW(IJ))
  175  CONTINUE

C*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
C*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**

c$doacross local(ij)
       DO 200 IJ=ISTART,IJLES
	  Y1(IJ)=ABS( VR(IJ)-VS(IJ) )
	  Y2(IJ)=ZR(IJ)*ZS(IJ)
	  Y3(IJ)=5./Y2(IJ)
	  Y4(IJ)=.08*Y3(IJ)*Y3(IJ)
	  Y5(IJ)=.05*Y3(IJ)*Y4(IJ)

	PRACS(IJ)=R7RF*Y1(IJ)*(Y3(IJ)/ZS(IJ)**5+Y4(IJ)/ZS(IJ)**3+Y5(IJ)
     1                        /ZS(IJ))
	QRACS(IJ)=MIN(D2T*PRACS(IJ), QS(IJ))
	PSACR(IJ)=R8RF*Y1(IJ)*(Y3(IJ)/ZR(IJ)**5+Y4(IJ)/ZR(IJ)**3+Y5(IJ)
     1                        /ZR(IJ))
	QSACR(IJ)=PSACR(IJ)

	 IF (TAIR(IJ) .GE. T0) THEN
	  PGAUT(IJ)=0.0
	  PRACS(IJ)=0.0
	  PSACR(IJ)=0.0
	 ELSE
	  QSACR(IJ)=0.0
	  QRACS(IJ)=0.0
	 ENDIF
C*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
C* 18 * PGFR : FREEZING OF QR TO QG                               **18**
	PGFR(IJ)=0.0
	PGAUT(IJ)=0.0
	 IF (TAIR(IJ) .LT. T0) THEN
	   Y1(IJ)=EXP(RN18A*(T0-TAIR(IJ)))
	  PGFR(IJ)=MAX(R18R*(Y1(IJ)-1.)/ZR(IJ)**7., 0.0)
	 ENDIF
  200 CONTINUE

C********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
C********   HANDLING THE NEGATIVE SNOW (QS)          *******************

c$doacross local(ij)
      DO 225 IJ=ISTART,IJLES
	  Y1(IJ)=QR(IJ)/D2T
	 PIACR(IJ)=MIN(Y1(IJ), PIACR(IJ))
	 DGACR(IJ)=MIN(Y1(IJ), DGACR(IJ))
	 PSACR(IJ)=MIN(Y1(IJ), PSACR(IJ))
	 PGFR(IJ)= MIN(Y1(IJ), PGFR(IJ))
	 Y1(IJ)=(PIACR(IJ)+DGACR(IJ)+PSACR(IJ)+PGFR(IJ))*D2T
	QR(IJ)=QR(IJ)+PR(IJ)+qracs(ij)-Y1(IJ)
	if (qr(ij) .lt. 0.0) then
	  y2(ij)=1.
	   if (y1(ij) .ne. 0.0) y2(ij)=qr(ij)/y1(ij)+1.
	  piacr(ij)=piacr(ij)*y2(ij)
	  dgacr(ij)=dgacr(ij)*y2(ij)
	  pgfr(ij)=pgfr(ij)*y2(ij)
	  psacr(ij)=psacr(ij)*y2(ij)
	  qr(ij)=0.0
	endif
	DLT2(IJ)=1.
	 IF (QR(IJ) .GT. 1.E-4) DLT2(IJ)=0.
	 IF (QS(IJ) .LE. 1.E-4) DLT2(IJ)=1.
	 IF (TAIR(IJ) .GE. T0) DLT2(IJ)=0.
	  Y1(IJ)=QS(IJ)/D2T
	 PGACS(IJ)=MIN(Y1(IJ), PGACS(IJ))
	 DGACS(IJ)=MIN(Y1(IJ), DGACS(IJ))
	 WGACS(IJ)=MIN(Y1(IJ), WGACS(IJ))
	 PGAUT(IJ)=MIN(Y1(IJ), PGAUT(IJ))
	 PRACS(IJ)=MIN(Y1(IJ), PRACS(IJ))
	 PWACS(IJ)=MIN(Y1(IJ), PWACS(IJ))
	PRN(IJ)=D2T*((1.-DLT3(IJ))*PIACR(IJ)+DGACR(IJ)+PGFR(IJ)
     1               +(1.-DLT2(IJ))*PSACR(IJ))
	PS(IJ)=PS(IJ)+D2T*(DLT3(IJ)*PIACR(IJ)+DLT2(IJ)*PSACR(IJ))
	 PRACS(IJ)=(1.-DLT2(IJ))*PRACS(IJ)
	 PWACS(IJ)=(1.-DLT4(IJ))*PWACS(IJ)

      PSN(IJ)=D2T*(PGACS(IJ)+DGACS(IJ)+WGACS(IJ)+PGAUT(IJ)+PRACS(IJ)
     1              +PWACS(IJ))

	QS(IJ)=QS(IJ)+PS(IJ)-qracs(ij)-PSN(IJ)
	 if (qs(ij) .lt. 0.0) then
	   y2(ij)=1.
	    if (psn(ij) .ne. 0.) y2(ij)=qs(ij)/psn(ij)+1.
	   pgacs(ij)=pgacs(ij)*y2(ij)
	   dgacs(ij)=dgacs(ij)*y2(ij)
	   wgacs(ij)=wgacs(ij)*y2(ij)
	   pgaut(ij)=pgaut(ij)*y2(ij)
	   pracs(ij)=pracs(ij)*y2(ij)
	   pwacs(ij)=pwacs(ij)*y2(ij)
	   qs(ij)=0.0
	 endif
      PSN(IJ)=D2T*(PGACS(IJ)+DGACS(IJ)+WGACS(IJ)+PGAUT(IJ)+PRACS(IJ)
     1             +PWACS(IJ))
       QG(IJ)=QG(IJ)+PG(IJ)+PRN(IJ)+PSN(IJ)
	Y1(IJ)=D2T*(PSACW(IJ)+PSFW(IJ)+DGACW(IJ)+PIACR(IJ)+DGACR(IJ)
     1             +PSACR(IJ)+PGFR(IJ))-QRACS(IJ)
       PT(IJ)=PT(IJ)+AFCP*Y1(IJ)
  225  CONTINUE
C* 11 * PSMLT : MELTING OF QS                                     **11**
C* 19 * PGMLT : MELTING OF QG TO QR                               **19**

c$doacross local(ij)
      DO 250 IJ=ISTART,IJLES
	  PSMLT(IJ)=0.0
	  PGMLT(IJ)=0.0
	TAIR(IJ)=(PT(IJ)+TB0)*PI0
	IF (TAIR(IJ) .GE. T0) THEN
	  TAIRC(IJ)=TAIR(IJ)-T0
	  DD(IJ)=R11T*TAIRC(IJ)*(R101R/ZS(IJ)**2+R102RF/ZS(IJ)**BSH5)
	 PSMLT(IJ)=MIN(QS(IJ),MAX(DD(IJ),0.0))
	   Y2(IJ)=R191R/ZG(IJ)**2+R192RF/ZG(IJ)**BGH5
	  DD1(IJ)=TAIRC(IJ)*(R19T*Y2(IJ)+R19AT*(QGACW(IJ)+QGACR(IJ)))
	 PGMLT(IJ)=MIN(QG(IJ),MAX(DD1(IJ),0.0))
	 PT(IJ)=PT(IJ)-AFCP*(PSMLT(IJ)+PGMLT(IJ))
	 QR(IJ)=QR(IJ)+PSMLT(IJ)+PGMLT(IJ)
	 QS(IJ)=QS(IJ)-PSMLT(IJ)
	 QG(IJ)=QG(IJ)-PGMLT(IJ)
	ENDIF
C* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
C* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0 < T <= T00)     **25**
C* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
	IF (QC(IJ).LE.CMIN) QC(IJ)=0.0
	IF (QI(IJ).LE.CMIN) QI(IJ)=0.0
	  TAIR(IJ)=(PT(IJ)+TB0)*PI0
	 IF(TAIR(IJ).LE.T00) THEN
	  PIHOM(IJ)=QC(IJ)
	 ELSE
	  PIHOM(IJ)=0.0
	 ENDIF
	 IF(TAIR(IJ).GE.T0) THEN
	  PIMLT(IJ)=QI(IJ)
	 ELSE
	  PIMLT(IJ)=0.0
	 ENDIF
	 PIDW(IJ)=0.0
	 IF (TAIR(IJ).LT.T0 .AND. TAIR(IJ).GT.T00) THEN
	   TAIRC(IJ)=TAIR(IJ)-T0
	    Y1(IJ)=MAX( MIN(TAIRC(IJ), -1.), -31.)
	    IT(IJ)=INT(ABS(Y1(IJ)))
C
	   Y2(IJ)=AA1(IT(IJ))
	   Y3(IJ)=AA2(IT(IJ))
	   Y4(IJ)=EXP(ABS(.5*TAIRC(IJ)))
	   DD(IJ)=(R00*QI(IJ)/(R25A*Y4(IJ)))**Y3(IJ)
	  PIDW(IJ)=MIN(R25RT*Y2(IJ)*Y4(IJ)*DD(IJ), QC(IJ))
	 ENDIF
	  Y1(IJ)=PIHOM(IJ)-PIMLT(IJ)+PIDW(IJ)
	PT(IJ)=PT(IJ)+AFCP*Y1(IJ)
	QC(IJ)=QC(IJ)-Y1(IJ)
	QI(IJ)=QI(IJ)+Y1(IJ)
C* 31 * PINT  : INITIATION OF QI                                  **31**
C* 32 * PIDEP : DEPOSITION OF QI                                  **32**
	PINT(IJ)=0.0
	 TAIR(IJ)=(PT(IJ)+TB0)*PI0
	 IF (TAIR(IJ) .LT. T0) THEN
	 if (qi(ij) .le. cmin2) qi(ij)=0.
	   TAIRC(IJ)=TAIR(IJ)-T0
	   DD(IJ)=R31R*EXP(BETA*TAIRC(IJ))
	    RTAIR(IJ)=1./(TAIR(IJ)-C76)
	     Y2(IJ)=EXP(C218-C580*RTAIR(IJ))
	    QSI(IJ)=RP0*Y2(IJ)
	    ESI(IJ)=C610*Y2(IJ)
	    SSI(IJ)=(QV(IJ)+QB0)/QSI(IJ)-1.
	   DM(IJ)=MAX( (QV(IJ)+QB0-QSI(IJ)), 0.)
	    RSUB1(IJ)=CS580*QSI(IJ)*RTAIR(IJ)*RTAIR(IJ)
	  DEP(IJ)=DM(IJ)/(1.+RSUB1(IJ))
	  PINT(IJ)=MAX(MIN(DD(IJ), DM(IJ)), 0.)
	     Y1(IJ)=1./TAIR(IJ)
	     Y2(IJ)=EXP(BETAH*TAIRC(IJ))
	     Y3(IJ)=SQRT(QI(IJ))
	   DD(IJ)=Y1(IJ)*(RN10A*Y1(IJ)-RN10B)+RN10C*TAIR(IJ)/ESI(IJ)
	  PIDEP(IJ)=MAX(R32RT*SSI(IJ)*Y2(IJ)*Y3(IJ)/DD(IJ), 0.)
	  PINT(IJ)=PINT(IJ)+PIDEP(IJ)
	   PINT(IJ)=MIN(PINT(IJ),DEP(IJ))
cc             if (pint(ij) .le. cmin2) pint(ij)=0.
	  PT(IJ)=PT(IJ)+ASCP*PINT(IJ)
	  QV(IJ)=QV(IJ)-PINT(IJ)
	  QI(IJ)=QI(IJ)+PINT(IJ)
	 ENDIF
  250  CONTINUE
C*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
        IF (NEW_ICE_SAT .EQ. 0) THEN ! cccshie by tao 5/3/01, shie 11/16/01 3d
c$doacross local(ij)
       DO 275 IJ=ISTART,IJLES
	 TAIR(IJ)=(PT(IJ)+TB0)*PI0
	CND(IJ)=RT0*(TAIR(IJ)-T00)
	DEP(IJ)=RT0*(T0-TAIR(IJ))
	  Y1(IJ)=1./(TAIR(IJ)-C358)
	  Y2(IJ)=1./(TAIR(IJ)-C76)
	 QSW(IJ)=RP0*EXP(C172-C409*Y1(IJ))
	 QSI(IJ)=RP0*EXP(C218-C580*Y2(IJ))
	  DD(IJ)=CP409*Y1(IJ)*Y1(IJ)
	  DD1(IJ)=CP580*Y2(IJ)*Y2(IJ)
	 IF (QC(IJ).LE.CMIN) QC(IJ)=CMIN
	 IF (QI(IJ).LE.CMIN) QI(IJ)=CMIN
	 IF (TAIR(IJ).GE.T0) THEN
	  DEP(IJ)=0.0
	  CND(IJ)=1.
	  QI(IJ)=0.0
	 ENDIF
	 IF (TAIR(IJ).LT.T00) THEN
	  CND(IJ)=0.0
	  DEP(IJ)=1.
	  QC(IJ)=0.0
	 ENDIF
	  Y5(IJ)=AVCP*CND(IJ)+ASCP*DEP(IJ)
	   Y1(IJ)=QC(IJ)*QSW(IJ)/(QC(IJ)+QI(IJ))
	   Y2(IJ)=QI(IJ)*QSI(IJ)/(QC(IJ)+QI(IJ))
	  Y4(IJ)=DD(IJ)*Y1(IJ)+DD1(IJ)*Y2(IJ)
	 QVS(IJ)=Y1(IJ)+Y2(IJ)
	 RSUB1(IJ)=(QV(IJ)+QB0-QVS(IJ))/(1.+Y4(IJ)*Y5(IJ))
	CND(IJ)=CND(IJ)*RSUB1(IJ)
	DEP(IJ)=DEP(IJ)*RSUB1(IJ)
	 IF (QC(IJ).LE.CMIN) QC(IJ)=0.
	 IF (QI(IJ).LE.CMIN) QI(IJ)=0.
CC    ******   CONDENSATION OR EVAPORATION OF QC  ******
	 CND(IJ)=MAX(-QC(IJ),CND(IJ))
CC    ******   DEPOSITION OR SUBLIMATION OF QI    ******
	 DEP(IJ)=MAX(-QI(IJ),DEP(IJ))
	PT(IJ)=PT(IJ)+AVCP*CND(IJ)+ASCP*DEP(IJ)
	QV(IJ)=QV(IJ)-CND(IJ)-DEP(IJ)
	QC(IJ)=QC(IJ)+CND(IJ)
	QI(IJ)=QI(IJ)+DEP(IJ)
  275  CONTINUE
       ENDIF ! cccshie by tao 5/3/01, shie 11/16/01 3d
C* 10 * PSDEP : DEPOSITION OF QS                                  **10**
C* 20 * PGDEP : DEPOSITION OF QG                                  **20**

cccshie 11/16/01 !  call tao et al (1989) saturation technique twice
        IF (NEW_ICE_SAT .EQ. 1) THEN
c$doacross local(ij)
           DO IJ=ISTART,IJLES
              TAIR(IJ)=(PT(IJ)+TB0)*PI0
c             CND1(IJ)=RT0*(TAIR(IJ)-T00)
c             DEP1(IJ)=RT0*(T0-TAIR(IJ))
              CND(IJ)=RT0*(TAIR(IJ)-T00)  ! cccshie by tao 5/3/01
              DEP(IJ)=RT0*(T0-TAIR(IJ))   ! cccshie by tao 5/3/01
                Y1(IJ)=1./(TAIR(IJ)-C358)
                Y2(IJ)=1./(TAIR(IJ)-C76)
              QSW(IJ)=RP0*EXP(C172-C409*Y1(IJ))
              QSI(IJ)=RP0*EXP(C218-C580*Y2(IJ))
              DD(IJ)=CP409*Y1(IJ)*Y1(IJ)
              DD1(IJ)=CP580*Y2(IJ)*Y2(IJ)
              Y5(IJ)=AVCP*CND(IJ)+ASCP*DEP(IJ) ! cccshie by tao 5/3/01
              Y1(IJ)=RT0*(TAIR(IJ)-T00)*QSW(IJ) ! cccshie by tao 5/3/01
              Y2(IJ)=RT0*(T0-TAIR(IJ))*QSI(IJ)  ! cccshie by tao 5/3/01
c             IF (QC(IJ).LE.CMIN) QC(IJ)=CMIN  ! cccshie by tao 5/3/01
c             IF (QI(IJ).LE.CMIN) QI(IJ)=CMIN  ! cccshie by tao 5/3/01
              IF (TAIR(IJ).GE.T0) THEN
c                 DEP1(IJ)=0.0
c                 CND1(IJ)=1.
c                 QI(IJ)=0.0
                  DEP(IJ)=0.0 ! cccshie by tao 5/3/01
                  CND(IJ)=1.  ! cccshie by tao 5/3/01
                  Y2(IJ)=0.   ! cccshie by tao 5/3/01
                  Y1(IJ)=QSW(IJ) ! cccshie by tao 5/3/01
              ENDIF
              IF (TAIR(IJ).LT.T00) THEN
                 CND(IJ)=0.0  ! cccshie by tao 5/3/01
                 DEP(IJ)=1.   ! cccshie by tao 5/3/01
                 Y2(IJ)=QSI(IJ) ! cccshie by tao 5/3/01
                 Y1(IJ)=0.     ! cccshie by tao 5/3/01
c                CND1(IJ)=0.0
c                DEP1(IJ)=1.
c                QC(IJ)=0.0
              ENDIF
c             Y5(IJ)=AVCP*CND1(IJ)+ASCP*DEP1(IJ) ! cccshie by tao 5/3/01
c             Y1(IJ)=QC(IJ)*QSW(IJ)/(QC(IJ)+QI(IJ)) ! cccshie by tao 5/3/01
c             Y2(IJ)=QI(IJ)*QSI(IJ)/(QC(IJ)+QI(IJ)) ! cccshie by tao 5/3/01
              Y4(IJ)=DD(IJ)*Y1(IJ)+DD1(IJ)*Y2(IJ)
              QVS(IJ)=Y1(IJ)+Y2(IJ)
              RSUB1(IJ)=(QV(IJ)+QB0-QVS(IJ))/(1.+Y4(IJ)*Y5(IJ))
             CND(IJ)=CND(IJ)*RSUB1(IJ) ! cccshie by tao 5/3/01
             DEP(IJ)=DEP(IJ)*RSUB1(IJ) ! cccshie by tao 5/3/01
c            CND1(IJ)=CND1(IJ)*RSUB1(IJ)
c            DEP1(IJ)=DEP1(IJ)*RSUB1(IJ)
c             IF (QC(IJ).LE.CMIN) QC(IJ)=0. ! cccshie by tao 5/3/01
c             IF (QI(IJ).LE.CMIN) QI(IJ)=0. ! cccshie by tao 5/3/01
CC    ******   CONDENSATION OR EVAPORATION OF QC  ******
c            CND1(IJ)=MAX(-QC(IJ),CND1(IJ))
             CND(IJ)=MAX(-QC(IJ),CND(IJ)) ! cccshie by tao 5/3/01
CC    ******   DEPOSITION OR SUBLIMATION OF QI    ******
             DEP(IJ)=MAX(-QI(IJ),DEP(IJ)) ! cccshie by tao 5/3/01
             PT(IJ)=PT(IJ)+AVCP*CND(IJ)+ASCP*DEP(IJ) ! cccshie by tao 5/3/01
             QV(IJ)=QV(IJ)-CND(IJ)-DEP(IJ) ! cccshie by tao 5/3/01
             QC(IJ)=QC(IJ)+CND(IJ) ! cccshie by tao 5/3/01
             QI(IJ)=QI(IJ)+DEP(IJ) ! cccshie by tao 5/3/01
c            DEP1(IJ)=MAX(-QI(IJ),DEP1(IJ))
c            PT(IJ)=PT(IJ)+AVCP*CND1(IJ)+ASCP*DEP1(IJ)
c            QV(IJ)=QV(IJ)-CND1(IJ)-DEP1(IJ)
c            QC(IJ)=QC(IJ)+CND1(IJ)
c            QI(IJ)=QI(IJ)+DEP1(IJ)
           ENDDO
        ENDIF ! IF (NEW_ICE_SAT .EQ. 1)


c$doacross local(ij)
      DO 280 IJ=ISTART,IJLES
       PSDEP(IJ)=0.0
       PGDEP(IJ)=0.0
       TAIR(IJ)=(PT(IJ)+TB0)*PI0
	IF (TAIR(IJ) .LT. T0) THEN
	  IF(QC(IJ)+QI(IJ).GT.1.E-5) THEN
	   DLT1(IJ)=1.
	  ELSE
	   DLT1(IJ)=0.
	  ENDIF
	 RTAIR(IJ)=1./(TAIR(IJ)-C76)
	  Y2(IJ)=EXP(C218-C580*RTAIR(IJ))
	 QSI(IJ)=RP0*Y2(IJ)
	 ESI(IJ)=C610*Y2(IJ)
	 SSI(IJ)=DLT1(IJ)*((QV(IJ)+QB0)/QSI(IJ)-1.)
	  DM(IJ)=QV(IJ)+QB0-QSI(IJ)
	  RSUB1(IJ)=CS580*QSI(IJ)*RTAIR(IJ)*RTAIR(IJ)
	  DD1(IJ)=MAX(DM(IJ)/(1.+RSUB1(IJ)),0.0)
	   Y3(IJ)=1./TAIR(IJ)
	  DD(IJ)=Y3(IJ)*(RN10A*Y3(IJ)-RN10B)+RN10C*TAIR(IJ)/ESI(IJ)
	   Y4(IJ)=R10T*SSI(IJ)*(R101R/ZS(IJ)**2+R102RF/ZS(IJ)**BSH5)
     1                         /DD(IJ)
	 PSDEP(IJ)=MAX(Y4(IJ), 0.0)
	  DD(IJ)=Y3(IJ)*(RN20A*Y3(IJ)-RN20B)+RN10C*TAIR(IJ)/ESI(IJ)
	  Y2(IJ)=R191R/ZG(IJ)**2+R192RF/ZG(IJ)**BGH5
	 PGDEP(IJ)=MAX(R20T*SSI(IJ)*Y2(IJ)/DD(IJ),0.0)
C     ******************************************************************
	   Y1(IJ)=MIN(PSDEP(IJ)+PGDEP(IJ),DD1(IJ))
	  PGDEP(IJ)=Y1(IJ)-PSDEP(IJ)
	 PT(IJ)=PT(IJ)+ASCP*Y1(IJ)
	 QV(IJ)=QV(IJ)-Y1(IJ)
	 QS(IJ)=QS(IJ)+PSDEP(IJ)
	 QG(IJ)=QG(IJ)+PGDEP(IJ)
	ENDIF
C* 23 * ERN : EVAPORATION OF QR                                   **23**
	ERN(IJ)=0.0
	IF (QR(IJ) .GT. 0.0) THEN
	 TAIR(IJ)=(PT(IJ)+TB0)*PI0
	  RTAIR(IJ)=1./(TAIR(IJ)-C358)
	   Y2(IJ)=EXP( C172-C409*RTAIR(IJ) )
	  ESW(IJ)=C610*Y2(IJ)
	  QSW(IJ)=RP0*Y2(IJ)
	  SSW(IJ)=(QV(IJ)+QB0)/QSW(IJ)-1.
	  DM(IJ)=QV(IJ)+QB0-QSW(IJ)
	   RSUB1(IJ)=CV409*QSW(IJ)*RTAIR(IJ)*RTAIR(IJ)
	  DD1(IJ)=MAX(-DM(IJ)/(1.+RSUB1(IJ)),0.0)
C           Y1(IJ)=R00*QRN(IJ,K)
C         ERN(IJ)=(((1.6+124.9*Y1(IJ)**.2046)*Y1(IJ)**.525)
C    1          /(2.55E6/(P00*QSW(IJ))+5.4E5))*(-DM(IJ)/(R00*QSW(IJ)))
C    2          *D2T
	    Y3(IJ)=1./TAIR(IJ)
	   DD(IJ)=Y3(IJ)*(RN30A*Y3(IJ)-RN10B)+RN10C*TAIR(IJ)/ESW(IJ)
	  y1(IJ)=-R23T*SSW(IJ)*(R231R/ZR(IJ)**2+R232RF/ZR(IJ)**3)
     1                          /DD(IJ)
	  ERN(IJ)=MIN(DD1(IJ),QR(IJ),MAX(y1(IJ),0.0))
	  PT(IJ)=PT(IJ)-AVCP*ERN(IJ)
	  QV(IJ)=QV(IJ)+ERN(IJ)
	  QR(IJ)=QR(IJ)-ERN(IJ)
	 ENDIF
  280  CONTINUE
C* 30 * PMLTG : EVAPORATION OF MELTING QG                         **30**
C* 33 * PMLTS : EVAPORATION OF MELTING QS                         **33**

c$doacross local(ij)
      DO 300 IJ=ISTART,IJLES
	PMLTS(IJ)=0.0
	PMLTG(IJ)=0.0
	 TAIR(IJ)=(PT(IJ)+TB0)*PI0
	IF (TAIR(IJ) .GE. T0) THEN
C           RTAIR(IJ)=1./(TAIR(IJ)-C358)
	   RTAIR(IJ)=1./(T0-C358)
	    Y2(IJ)=EXP( C172-C409*RTAIR(IJ) )
	   ESW(IJ)=C610*Y2(IJ)
	   QSW(IJ)=RP0*Y2(IJ)
	   SSW(IJ)=1.-(QV(IJ)+QB0)/QSW(IJ)
	   DM(IJ)=QSW(IJ)-QV(IJ)-QB0
	   RSUB1(IJ)=CV409*QSW(IJ)*RTAIR(IJ)*RTAIR(IJ)
	  DD1(IJ)=MAX(DM(IJ)/(1.+RSUB1(IJ)),0.0)
	    Y3(IJ)=1./TAIR(IJ)
	   DD(IJ)=Y3(IJ)*(RN30A*Y3(IJ)-RN10B)+RN10C*TAIR(IJ)/ESW(IJ)
	   Y1(IJ)=R30T*SSW(IJ)*(R191R/ZG(IJ)**2+R192RF/ZG(IJ)**BGH5)
     1                          /DD(IJ)
	  PMLTG(IJ)=MIN(QG(IJ),MAX(Y1(IJ),0.0))
	   Y1(IJ)=R33T*SSW(IJ)*(R331R/ZS(IJ)**2+R332RF/ZS(IJ)**BSH5)
     1                          /DD(IJ)
	  PMLTS(IJ)=MIN(QS(IJ),MAX(Y1(IJ),0.0))
	   Y1(IJ)=MIN(PMLTG(IJ)+PMLTS(IJ),DD1(IJ))
	  PMLTG(IJ)=Y1(IJ)-PMLTS(IJ)
	  PT(IJ)=PT(IJ)-ASCP*Y1(IJ)
	  QV(IJ)=QV(IJ)+Y1(IJ)
	  QS(IJ)=QS(IJ)-PMLTS(IJ)
	  QG(IJ)=QG(IJ)-PMLTG(IJ)
	ENDIF
C       IF (QV(IJ)+QB0 .LE. 0.) QV(IJ)=-QB0
	 IF(QC(IJ).LT.CMIN) QC(IJ)=0.0
	 IF(QR(IJ).LT.CMIN) QR(IJ)=0.0
	 IF(QI(IJ).LT.CMIN) QI(IJ)=0.0
	 IF(QS(IJ).LT.CMIN) QS(IJ)=0.0
	 IF(QG(IJ).LT.CMIN) QG(IJ)=0.0

       DPT(IJ,K)=PT(IJ)
       DQV(IJ,K)=QV(IJ)
       QCL(IJ,K)=QC(IJ)
       QRN(IJ,K)=QR(IJ)
       QCI(IJ,K)=QI(IJ)
       QCS(IJ,K)=QS(IJ)
       QCG(IJ,K)=QG(IJ)

c      q1_hyd(ij,k)=q1_hyd(ij,k)+pt(ij)-tttbud(ij)
c      q2_hyd(ij,k)=q2_hyd(ij,k)+qv(ij)-qqqbud(ij)

c      if(isec.eq.isec/ibudsec*ibudsec) then
c         q1a_hyd(ij,k) = q1_hyd(ij,k) / rbud
c         q2a_hyd(ij,k) = q2_hyd(ij,k) / rbud
c         q1_hyd(ij,k)=0.
c         q2_hyd(ij,k)=0.
c      endif

  300  CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   sddd=0.
	   ssss=0.
	   shhh=0.
	   sccc=0.
	   smmm=0.
	   sfff=0.

       DO 305 IJ=ISTART,IJLES
	 DD(IJ)=MAX(-CND(IJ), 0.)
	 CND(IJ)=MAX(CND(IJ), 0.)
	 DD1(IJ)=MAX(-DEP(IJ), 0.)
	 DEP(IJ)=MAX(DEP(IJ), 0.)
	 SCC=SCC+CND(IJ)*ASSS(IJ)
	 SEE=SEE+(DD(IJ)+ERN(IJ))*ASSS(IJ)
c
	 sddd=sddd+(dep(ij)+pint(ij)+psdep(ij)+pgdep(ij))*asss(ij)
	 ssss=ssss+dd1(ij)*asss(ij)
	 shhh=shhh+rsw(ij,k)*dt*asss(ij)
	 sccc=sccc+rlw(ij,k)*dt*asss(ij)
	 smmm=smmm+(psmlt(ij)+pgmlt(ij)+pimlt(ij))*asss(ij)
	 sfff=sfff+asss(ij)*dt*(psacw(ij)+piacr(ij)+psfw(ij)+pgfr(ij)
     1        +dgacw(ij)+dgacr(ij)+psacr(ij))
     2        -qracs(ij)+pihom(ij)+pidw(ij)
c
  305   CONTINUE
c
	s_dep(k)=s_dep(k)+sddd
	s_sub(k)=s_sub(k)+ssss
	s_qrs(k)=s_qrs(k)+shhh
	s_qrl(k)=s_qrl(k)+sccc
	s_mel(k)=s_mel(k)+smmm
	s_frz(k)=s_frz(k)+sfff
c
	SC(K)=SCC+SC(K)
	SE(K)=SEE+SE(K)
CC    ***   STATISTICS FOR CONVECTIVE AND ANVIL REGIMES   ***********

       IF (ID .EQ. 1) THEN
	 RDTS=1./D2T

c$doacross local(ij)
	DO 310 IJ=ISTART,IJLES
	   CND(IJ)=CND(IJ)*RDTS
	   DD(IJ)=DD(IJ)*RDTS
	   PINT(IJ)=PINT(IJ)*RDTS
	   PIDW(IJ)=PIDW(IJ)*RDTS
	   PIMLT(IJ)=PIMLT(IJ)*RDTS
	   PIHOM(IJ)=PIHOM(IJ)*RDTS
	   PSMLT(IJ)=PSMLT(IJ)*RDTS
	   PGMLT(IJ)=PGMLT(IJ)*RDTS
	   PSDEP(IJ)=PSDEP(IJ)*RDTS
	   PGDEP(IJ)=PGDEP(IJ)*RDTS
	   PMLTG(IJ)=PMLTG(IJ)*RDTS
	   PMLTS(IJ)=PMLTS(IJ)*RDTS
	   DD1(IJ)=DD1(IJ)*RDTS
	   DEP(IJ)=DEP(IJ)*RDTS
	   ERN(IJ)=ERN(IJ)*RDTS
	   QRACS(IJ)=QRACS(IJ)*RDTS
	  DDA(IJ)=RSW(IJ,K)
	  DDB(IJ)=RLW(IJ,K)
	 Y1(IJ)=QC(IJ)+QR(IJ)+QI(IJ)+QS(IJ)+QG(IJ)
c
	 dm(ij)=a0*(rho1(k)*ww1(ij,k)+rho1(k+1)*ww1(ij,k+1)+
     1              y0(ij)*(rho1(k)*wb(k)+rho1(k+1)*wb(k+1)))
	 rq(ij)=.005*(rho1(k)*(ww1(ij,k)+wb(k))+
     1               rho1(k+1)*(ww1(ij,k+1)+wb(k+1)))/r00
  310    continue


C       DO 1050 KC=1,7
	KC=4
	 DO MT=1,4
c$doacross local(ij)
	 DO IJ=ISTART,IJLES
	  IBZ(IJ,MT)=0
	   IF(ICS5(IJ,MT).EQ.1) IBZ(IJ,MT)=1
	 ENDDO
	 ENDDO

c$doacross local(ij)
	DO 315 IJ=ISTART,IJLES
	    IBZ(IJ,1)=1
  315    CONTINUE

	 DO MT=1,4

	 DO 330 IJ=ISTART,IJLES
	   IF(KC.EQ.4) GO TO 330
	   IF(KC.LE.3) GO TO 36
	    IF (RQ(IJ).GT.RBY(KC)) IBZ(IJ,MT)=0
	    GO TO 330
   36       IF (RQ(IJ).LT.RBY(KC)) IBZ(IJ,MT)=0
  330    CONTINUE
	 ENDDO

c$doacross local(mt,ij,sww,scc,see,a1,a2,a3,a4,a5,a6,a7,a8,a9,
c$&    a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,
c$&    a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,
c$&    a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50)
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
	  A45=0.
	  A46=0.
	  A47=0.
	  A48=0.
	  A49=0.
	  A50=0.

cc$doacross local(ij),reduction(sww,scc,see,a1,a2,a3,a4,a5,a6,a7,a8,a9,
cc$&    a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25)
          DO 30 IJ=ISTART,IJLES
	     ASSS(IJ)=CSTTT(IJ)
	     IF(IBZ(IJ,MT).EQ.1) THEN
	        SWW=SWW+DM(IJ)*ASSS(IJ)
	        SCC=SCC+CND(IJ)*ASSS(IJ)
	        SEE=SEE+(DD(IJ)+ERN(IJ))*ASSS(IJ)
	        a1=a1+(pihom(ij)+pidw(ij))*asss(ij)
	        A2=A2+PINT(IJ)*ASSS(IJ)
	        A3=A3+PGFR(IJ)*ASSS(IJ)
	        A4=A4+PSAUT(IJ)*ASSS(IJ)
	        A5=A5+PSACI(IJ)*ASSS(IJ)
	        A6=A6+PSACW(IJ)*ASSS(IJ)
	        A7=A7+PRACI(IJ)*ASSS(IJ)
	        A8=A8+PIACR(IJ)*ASSS(IJ)
	        A9=A9+PRAUT(IJ)*ASSS(IJ)
	        A10=A10+PRACW(IJ)*ASSS(IJ)
	        A11=A11+PSFW(IJ)*ASSS(IJ)
	        A12=A12+PSFI(IJ)*ASSS(IJ)
	        A13=A13+(PGACS(IJ)+DGACS(IJ))*ASSS(IJ)
	        A14=A14+DGACW(IJ)*ASSS(IJ)
	        A15=A15+DGACI(IJ)*ASSS(IJ)
	        A16=A16+DGACR(IJ)*ASSS(IJ)
	        A17=A17+PMLTG(IJ)*ASSS(IJ)
	        A18=A18+DEP(IJ)*ASSS(IJ)
	        A19=A19+PRACS(IJ)*ASSS(IJ)
	        A20=A20+PSACR(IJ)*ASSS(IJ)
	        A21=A21+PMLTS(IJ)*ASSS(IJ)
	        A22=A22+PSMLT(IJ)*ASSS(IJ)
	        A23=A23+PGMLT(IJ)*ASSS(IJ)
	        A24=A24+PSDEP(IJ)*ASSS(IJ)
	        A25=A25+PIMLT(IJ)*ASSS(IJ)
	     ENDIF
   30     CONTINUE

cc$doacross local(ij),reduction(a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,
cc$&    a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50)
	 DO 35 IJ=ISTART,IJLES
	  IF(IBZ(IJ,MT).EQ.1) THEN
	  A26=A26+PGDEP(IJ)*ASSS(IJ)
	  A27=A27+DD1(IJ)*ASSS(IJ)
	  A28=A28+PRACI(IJ)*DLT3(IJ)*ASSS(IJ)
	  A29=A29+PIACR(IJ)*DLT3(IJ)*ASSS(IJ)
	  A30=A30+PSACR(IJ)*DLT2(IJ)*ASSS(IJ)
	  A31=A31+QRACS(IJ)*ASSS(IJ)
	  A32=A32+PSACW(IJ)*DLT4(IJ)*ASSS(IJ)
	  A33=A33+QCL(IJ,K)*ASSS(IJ)
	  A34=A34+QRN(IJ,K)*ASSS(IJ)
	  A35=A35+QCI(IJ,K)*ASSS(IJ)
	  A36=A36+QCS(IJ,K)*ASSS(IJ)
	  A37=A37+QCG(IJ,K)*ASSS(IJ)
	  A38=A38+ERN(IJ)*ASSS(IJ)
	  A39=A39+WGACR(IJ)*ASSS(IJ)
	  A40=A40+QSACW(IJ)*ASSS(IJ)
	  A41=A41+DDA(IJ)*ASSS(IJ)
	  A42=A42+DDB(IJ)*ASSS(IJ)
	   A43=A43+(QV(IJ)+QA1(K)-QA(K))*ASSS(IJ)
	   A44=A44+(PT(IJ)+TA1(K)-TA(K))*ASSS(IJ)
	   A45=A45+ASSS(IJ)
	    A46=A46+Y1(IJ)*ASSS(IJ)
	  A47=A47+(PSACW(IJ)+PSFW(IJ)+DGACW(IJ)+PIACR(IJ)+DGACR(IJ)
     1        +PSACR(IJ)+PGFR(IJ)-QRACS(IJ)+PIHOM(IJ)-PIMLT(IJ)
     2        +PIDW(IJ))*ASSS(IJ)
	  A48=A48+(Y1(IJ)-QCL1(IJ,K)-QRN1(IJ,K)-QCI1(IJ,K)-QCS1(IJ,K)
     1            -QCG1(IJ,K))*ASSS(IJ)
	  A49=A49+(QV(IJ)-DQV1(IJ,K))*ASSS(IJ)
	  A50=A50+(PT(IJ)-DPT1(IJ,K))*ASSS(IJ)
	 ENDIF
   35   CONTINUE

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
	Q1T(K,MT,KC)=Q1T(K,MT,KC)+A47
	QC0(K,MT,KC)=A33
	QR0(K,MT,KC)=A34
	QI0(K,MT,KC)=A35
	QS0(K,MT,KC)=A36
	QG0(K,MT,KC)=A37
	QV0(K,MT,KC)=A43
	TT0(K,MT,KC)=A44
	SGPT(K,MT,KC)=A45
	TSQQ(K,MT,KC)=A46
	SQHDT(K,MT,KC)=SQHDT(K,MT,KC)+A48
	SQVDT(K,MT,KC)=SQVDT(K,MT,KC)+A49
	SQTDT(K,MT,KC)=SQTDT(K,MT,KC)+A50
	SQC0(K,MT,KC)=SQC0(K,MT,KC)+QC0(K,MT,KC)
	SQR0(K,MT,KC)=SQR0(K,MT,KC)+QR0(K,MT,KC)
	SQI0(K,MT,KC)=SQI0(K,MT,KC)+QI0(K,MT,KC)
	SQS0(K,MT,KC)=SQS0(K,MT,KC)+QS0(K,MT,KC)
	SQG0(K,MT,KC)=SQG0(K,MT,KC)+QG0(K,MT,KC)
	SQV0(K,MT,KC)=SQV0(K,MT,KC)+QV0(K,MT,KC)
	STT0(K,MT,KC)=STT0(K,MT,KC)+TT0(K,MT,KC)
	SSGPT(K,MT,KC)=SSGPT(K,MT,KC)+SGPT(K,MT,KC)
	TSQQ1(K,MT,KC)=TSQQ1(K,MT,KC)+TSQQ(K,MT,KC)
  350   CONTINUE
       ENDIF
CC    ********************************************
       IF(ID.EQ.1) THEN
C--------------------------------------------------------------------
C    CONDENSATION:  CND(IJ)
C    EVAPORATION:   DD(IJ)+ERN(IJ)
C    DEPOSITION:    DEP(IJ)+PSDEP(IJ)+PGDEP(IJ)+PINT(IJ)
C    SUBLIMATION:   DD1(IJ)+PMLTS(IJ)+PMLTG(IJ)
C    MELTING:       PSMLT(IJ)+PGMLT(IJ)+PIMLT(IJ)+QRACS(IJ)
C    FREEZING:      PIHOM(IJ)+PIDW(IJ)+PSACW(IJ)+PSFW(IJ)+DGACW(IJ)
C                   +PIACR(IJ)+DGACR(IJ)+PSACR(IJ)+PGFR(IJ)
C    MASS FLUX:     DM(IJ)
C    CLOUD WATER:   QC(IJ)
C    RAIN:          
C    CLOUD ICE
C    SNOW
C    HAIL/GRAUPEL:
C----------------------------------------------------------------------
c$doacross local(ij,a1,a2,a3,a11,a22,a33,zdry,a44,zwet)
       DO 42 IJ=ISTART,IJLES
	CND(IJ)=CND(IJ)*RFT*ASSS(IJ)
	ERN(IJ)=(ERN(IJ)+DD(IJ))*RFT*ASSS(IJ)
	Y1(IJ)=(DEP(IJ)+PSDEP(IJ)+PGDEP(IJ)+PINT(IJ))*RFT*ASSS(IJ)
	Y2(IJ)=(DD1(IJ)+PMLTS(IJ)+PMLTG(IJ))*RFT*ASSS(IJ)
	Y3(IJ)=(PSMLT(IJ)+PGMLT(IJ)+PIMLT(IJ)+QRACS(IJ))*RFT*ASSS(IJ)
	Y4(IJ)=(PIHOM(IJ)+PIDW(IJ)+PSACW(IJ)+PSFW(IJ)+DGACW(IJ)+
     1          PIACR(IJ)+DGACR(IJ)+PSACR(IJ)+PGFR(IJ)+(pihom(ij)+
     2          pidw(ij))*rdt)*RFT*ASSS(IJ)
	Y5(IJ)=DM(IJ)*ASSS(IJ)
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
	DBZ(IJ)=DBZ(IJ)*ASSS(IJ)
	QC(IJ)=QC(IJ)*ASSS(IJ)
	QR(IJ)=QR(IJ)*ASSS(IJ)
	QI(IJ)=QI(IJ)*ASSS(IJ)
	QS(IJ)=QS(IJ)*ASSS(IJ)
	QG(IJ)=QG(IJ)*ASSS(IJ)

   42  CONTINUE

	DO 44 IJ=ISTART,IJLES
	IF(RQ(IJ) .GE. 0.) SCU1(K)=SCU1(K)+CND(IJ)
	IF(RQ(IJ) .LT. 0.) SED1(K)=SED1(K)+ERN(IJ)
   44  CONTINUE

       DO 40 IJ=ISTART,IJLES
	 if(rq(ij) .gt. -0.5) then
            iww=min(rq(ij)+0.5, 15.)+1
	    IF (ICS5(IJ,2) .EQ. 1) THEN
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
	       SCNT(IWW,K)=SCNT(IWW,K)+ASSS(IJ)
	    ELSEIF (ICS5(IJ,3) .EQ. 1) THEN
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
	       SSCNT(IWW,K)=SSCNT(IWW,K)+ASSS(IJ)
	    ENDIF
	 else
            jww=min(-rq(ij)+0.5, 5.)
	    IF (ICS5(IJ,2) .EQ. 1) THEN
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
	       SNCNT(jww,K)=SNCNT(jww,K)+ASSS(IJ)
	    ELSEIF (ICS5(IJ,3) .EQ. 1) THEN
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
	       SSNCNT(jww,K)=SSNCNT(jww,K)+ASSS(IJ)
	    ENDIF
         endif
   40  CONTINUE

       ENDIF
 1000 CONTINUE
C     ****************************************************************
      IF (ID .EQ. 1) THEN
       LC=0
       LS=0
       DO 390 IJ=ISTART,IJLES
	LC=LC+ICS(IJ,2)
	LS=LS+ICS(IJ,3)
  390  CONTINUE
	IF (LC .EQ. 0) LC=1000000
	IF (LS .EQ. 0) LS=1000000
       DO 400 MT=1,3
	A1=RIJL2
	 IF (MT .EQ. 2) A1=1./FLOAT(LC)
	 IF (MT .EQ. 3) A1=1./FLOAT(LS)

       DO K=KT1,KT2
	B1(K)=0.0
	B2(K)=0.0
       DO IJ=ISTART,IJLES
	IF(ICS5(IJ,MT).EQ.1) THEN
	 B1(K)=B1(K)+DPT(IJ,K)
	 B2(K)=B2(K)+DQV(IJ,K)
	ENDIF
       ENDDO
       ENDDO

       DO 430 K=KT1,KT2
	 TB00(K,MT)=B1(K)*A1
	 QB00(K,MT)=B2(K)*A1
  430  CONTINUE
  400 CONTINUE
      ENDIF
      D2T=D22T
      CALL SATDT
      RETURN
      END    

      SUBROUTINE TERVRH (IRSG,RHO,FV)
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      PARAMETER (NXY=NX*NY,NM16=16*NM)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/B3CS/ AG,BG,AS,BS,AW,BW,BGH,BGQ,BSH,BSQ,BWH,BWQ
      COMMON/RTERV/ ZRC,ZGC,ZSC,VR0,VR1,VR2,VR3,VGC,VSC
      COMMON/B1CR/ QCL(NXY,NZ),QRN(NXY,NZ)
      COMMON/B1IG/ QCI(NXY,NZ),QCG(NXY,NZ)
      COMMON/B1S/ QCS(NXY,NZ)
      COMMON/BADV1/ WW1(NXY,NZ)
      COMMON/BBA/ Y1(NXY),VR(NXY),VS(NXY),VG(NXY),Z1(NM16)
      DIMENSION FV(1),RHO(1)
      save
C     ******************************************************************
c       IJLES=NX*(NY-1)-1
c       ISTART=NX+2
c
       ijles=max(nx*(ny-1), ny*(nx-1))-1
       istart=min(nx,ny)+2
c
       CMIN=1.E-15
c$doacross local(k,ij)
      DO K=1,KMAX
      DO IJ=1,NXY
	WW1(IJ,K)=0.
      ENDDO
      ENDDO
      IF(IRSG.NE.0) GO TO 1

c$doacross local(k,a1,km,ij)
      DO 100 K=2,KLES
	A1=.5*RHO(K)
	KM=K-1
       DO 10 IJ=ISTART,IJLES
	  Y1(IJ)=A1*(QRN(IJ,K)+QRN(IJ,KM))
	IF (Y1(IJ) .GT. CMIN) THEN
	  VS(IJ)=SQRT( Y1(IJ) )
	  VG(IJ)=SQRT( VS(IJ) )
	  VR(IJ)=VR0+VR1*VG(IJ)+VR2*VS(IJ)+VR3*VG(IJ)*VS(IJ)
	 WW1(IJ,K)=MAX(FV(K)*VR(IJ), 0.e0)
	ENDIF
   10  CONTINUE
  100 CONTINUE

      RETURN
    1 IF(IRSG.NE.1) GO TO 2
      DO 200 K=2,KLES
	 KM=K-1
	 A1=.5*RHO(K)
	VSCF=VSC*FV(K)
c$doacross local(ij)
       DO 20 IJ=ISTART,IJLES
	 Y1(IJ)=A1*(QCS(IJ,K)+QCS(IJ,KM))
	IF (Y1(IJ) .GT. CMIN) THEN
	 WW1(IJ,K)=MAX(VSCF*Y1(IJ)**BSQ, 0.e0)
	ENDIF
   20  CONTINUE
  200 CONTINUE
      RETURN
    2 DO 300 K=2,KLES
	KM=K-1
	 A1=.5*RHO(K)
	VGCR=VGC*FV(K)
c$doacross local(ij)
       DO 30 IJ=ISTART,IJLES
	  Y1(IJ)=A1*(QCG(IJ,K)+QCG(IJ,KM))
	IF (Y1(IJ) .GT. CMIN) THEN
	 WW1(IJ,K)=MAX(VGCR*Y1(IJ)**BGQ, 0.e0)
	ENDIF
   30  CONTINUE
  300 CONTINUE
      RETURN
      END
      SUBROUTINE SATICEW
C     (SOONG & OGURA, 1973)    COMPUTE WATER PHASE MICROPHYSICS AND 
C     SATURATION PROCESSES
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880)
      PARAMETER (NXY=NX*NY,NT2=2*NT)
      PARAMETER (NZ2=2*NZ,NZ3=3*NZ,NZ4=4*NZ,NM2=2*NM,NM11=11*NM)
      PARAMETER (NB=NX*NY*(NZ-29),NB1=NX*NY*(NZ-26))
      PARAMETER (NB3=NZ+NX,NB4=5*NT,NB5=NX*NY*(NZ-27))
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1   NRAN,KT1,KT2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1   PSFC(5)
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/B3CS/ AG,BG,AS,BS,AW,BW,BGH,BGQ,BSH,BSQ,BWH,BWQ
      COMMON/SIZE/ TNW,TNS,TNG,ROQS,ROQG,ROQR
      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC
      COMMON/BSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,BND1,RN2,BND2,
     1   RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,RN101,RN10A,RN11,RN11A,
     2   RN12,RN12A(31),RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN17,
     3   RN17A,RN17B,RN17C,RN18,RN18A,RN19,RN19A,RN19B,RN20,RN20A,RN20B,
     4   BND3,RN21,RN22,RN23,RN23A,RN23B,RN25,RN25A(31),RN30A,RN30B,
     5   RN30C,RN31,BETA,RN32
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/BSAT1/ PT(NXY),QV(NXY),QC(NXY),QR(NXY),QI(NXY),QS(NXY),
     1 QG(NXY),WGACR(NXY),DEP(NXY),DD(NXY),DD1(NXY),QVS(NXY),DM(NXY),
     2 RSUB1(NXY),COL(NXY),CND(NXY),RQ(NXY),ERN(NXY),SCV(NXY),TCA(NXY),
     3 DWV(NXY),ZR(NXY),VR(NXY),ZS(NXY),VS(NXY),ZG(NXY),VG(NXY),
     4 EGS(NXY),DBZ(NXY),DDD(NB)
      COMMON/BADV/ ESI(NXY),QSI(NXY),SSI(NXY),QSW(NXY),SSW(NXY),
     1 PIHOM(NXY),PIDW(NXY),PIMLT(NXY),PSAUT(NXY),PSACI(NXY),QSACW(NXY),
     2 PRACI(NXY),PIACR(NXY),PRAUT(NXY),PRACW(NXY),PSFW(NXY),PSFI(NXY),
     3 DGACS(NXY),DGACW(NXY),DGACI(NXY),DGACR(NXY),PGACS(NXY),
     4 WGACS(NXY),QGACW(NXY),WGACI(NXY),QGACR(NXY),PSACW(NXY),
     5 DDA(NB5)
      COMMON/BSAT/ PGAUT(NXY),PRACS(NXY),PSACR(NXY),QSACR(NXY),
     1 PGFR(NXY),PSMLT(NXY),PGMLT(NXY),PSDEP(NXY),PS(NXY),PSSUB(NXY),
     2 PGSUB(NXY),TAIR(NXY),TAIRC(NXY),PR(NXY),PG(NXY),PRN(NXY),
     3 PSN(NXY),DLT1(NXY),DLT2(NXY),DLT3(NXY),RTAIR(NXY),Y5(NXY),
     4 PINT(NXY),PIDEP(NXY),PGWET(NXY),ASSS(NXY),DDB(NB1)
      COMMON/B1TQ/ DPT(NXY,NZ),DQV(NXY,NZ)
      COMMON/B1CR/ QCL(NXY,NZ),QRN(NXY,NZ)
      COMMON/B2TQ/ DPT1(NXY,NZ),DQV1(NXY,NZ)
      COMMON/B2CR/ QCL1(NXY,NZ),QRN1(NXY,NZ)
      COMMON/B4WP/ WW1(NXY,NZ)
      COMMON/SLWAVE/ RSW(NX,NY,NZ),RLW(NX,NY,NZ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc      COMMON/BW/ TRAHX(NXY,NZ),TRAHY(NXY,NZ),TRAV(NXY,NZ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),ZX(NZ4),AM(NZ),ZQ(NZ3),WB(NZ),ZW(NZ2),RRHO(NZ),WBX(NB3)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      COMMON/BRH1/ SRRO(NZ),QRRO(NZ),SQC(NZ),SQR(NZ),SQI(NZ),SQS(NZ),
     1   SQG(NZ),STQC(NZ),STQR(NZ),STQI(NZ),STQS(NZ),STQG(NZ),TTTD(NB4)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      COMMON/BSTS/ THOM(NZ,4,7),TDW(NZ,4,7),TMLT(NZ,4,7),SAUT(NZ,4,7),
     1 SACI(NZ,4,7),SACW(NZ,4,7),RACI(NZ,4,7),TACR(NZ,4,7),RAUT(NZ,4,7),
     2 RACW(NZ,4,7),SFW(NZ,4,7),SFI(NZ,4,7),GACS(NZ,4,7),GACW(NZ,4,7),
     3 GACI(NZ,4,7),GACR(NZ,4,7),GWET(NZ,4,7),GAUT(NZ,4,7),RACS(NZ,4,7),
     4 SACR(NZ,4,7),GFR(NZ,4,7),SMLT(NZ,4,7),GMLT(NZ,4,7),SDEP(NZ,4,7),
     5 SSUB(NZ,4,7),GSUB(NZ,4,7),PERN(NZ,4,7),D3RI(NZ,4,7),D3IR(NZ,4,7),
     6 D2SR(NZ,4,7),D2RS(NZ,4,7),GDRY(NZ,4,7),COC(NZ,4,7),COE(NZ,4,7),
     7 SMF0(NZ,4,7),QC0(NZ,4,7),QR0(NZ,4,7),QI0(NZ,4,7),QS0(NZ,4,7),
     8 QG0(NZ,4,7),SQC0(NZ,4,7),SQR0(NZ,4,7),SQI0(NZ,4,7),SQS0(NZ,4,7),
     9 SQG0(NZ,4,7),ERNS(NZ,4,7),WGRS(NZ,4,7),QSWS(NZ,4,7),TB00(NZ,4),
     1 QB00(NZ,4)
      COMMON/BSTS1/ TUT1(NZ,4,7),TUT2(NZ,4,7),TVT1(NZ,4,7),TVT2(NZ,4,7),
     1 TSTF(NZ,4,7),TSTF1(NZ,4,7),TSTF2(NZ,4,7),TSQF(NZ,4,7),
     2 QQQ(NZ,4,7),TSQF1(NZ,4,7),TSQF2(NZ,4,7),TSQQ(NZ,4,7),
     3 TSQQ1(NZ,4,7)
      COMMON/BSTS3/ QV0(NZ,4,7),TT0(NZ,4,7),SQV0(NZ,4,7),STT0(NZ,4,7),
     1 SGPT(NZ,4,7),SSGPT(NZ,4,7),SNQHD(NZ,4,7),SNQVD(NZ,4,7),
     2 Q1T(NZ,4,7),SNHDH(NZ,4,7),SQHDT(NZ,4,7),SQVDT(NZ,4,7)
      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7),SQTDT(NZ,4,7),SQHL(NZ,4,7)

      common/bsts40/ fcld(nz,4,7)

      COMMON/BCS/ S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     1 S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     2 S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     3 SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     4 SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     5 SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCSS/ SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     1 SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     2 SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     3 SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     4 SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     5 SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)
      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NXY),ICS5(NXY,4),
     1  IBZ(NXY,4)
      COMMON/BI/ IT(NXY),ICS(NXY,4)
      COMMON/RSTAT/ CSTTT(NXY),CSTT(NXY)
C
      COMMON/BBA/ Y1(NXY),Y2(NXY),Y3(NXY),Y4(NXY),Y7(NM2),B0(NM),B1(NM),
     1   B2(NM),Y6(NM11)
C
c      DIMENSION FV(NZ)
      REAL RBY(7)
      DATA RBY/2.,1.,0.,0.,0.,-.5,-1./
c      save d22t
      save
CC    ***   TWO CLASSES OF WATER-PHASE   *******************************
       D22T=D2T
      IF(IJKADV .EQ. 1) THEN
	 D2T=DT
      ELSE
	 D2T=D2T
      ENDIF
c       IJLES=NX*(NY-1)-1
c       ISTART=NX+2
c
       ijles=max(nx*(ny-1), ny*(nx-1))-1
       istart=min(nx,ny)+2
c
       CMIN=1.E-20
       CMIN2=1.E-20
	UCOR=3071.29/TNW**.75                                            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO IJ=1,NXY
       IT(IJ)=1
      ENDDO
      F2=RD1*D2T
      F3=RD2*D2T
      FT=DT/D2T
      RFT=RIJL2*FT
      A0=.5*ISTATMIN*RIJL2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 1000 K=KT1,KT2
c        if (ijkadv .eq. 1) then
c          tb0=ta(k)
c          qb0=qa(k)
c        else
	  tb0=ta1(k)
	  qb0=qa1(k)
c        endif
       P00=P0(K)
       RP0=3.799052E3/P0(K)
       PI0=PI(K)
       R00=RHO(K)
       F00=F0(K)
CCCCCC
       DO 250 IJ=ISTART,IJLES
	PT(IJ)=DPT(IJ,K)
	QV(IJ)=DQV(IJ,K)
	QC(IJ)=QCL(IJ,K)
	QR(IJ)=QRN(IJ,K)
C        IF (QV(IJ)+QB0 .LE. 0.) QV(IJ)=-QB0
	IF (QC(IJ) .LE. CMIN) QC(IJ)=0.
	IF (QR(IJ) .LE. CMIN) QR(IJ)=0.
       RTAIR(IJ)=1./((PT(IJ)+TB0)*PI0-C358)
       DD(IJ)=C172-C409*RTAIR(IJ)
       QVS(IJ)=RP0*EXP(DD(IJ))
       DM(IJ)=QV(IJ)+QB0-QVS(IJ)
       RSUB1(IJ)=QVS(IJ)*F5*RTAIR(IJ)*RTAIR(IJ)
C     ******   AUTOCONVERSION AND COLLECTION   ******
       PRAUT(IJ)=MAX(F2*(QC(IJ)-BOUND), 0.0)
       PRACW(IJ)=F3*QC(IJ)*QR(IJ)**.875
C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
	   Y1(IJ)=QC(IJ)/D2T
	  PRAUT(IJ)=MIN(Y1(IJ), PRAUT(IJ))
	  PRACW(IJ)=MIN(Y1(IJ), PRACW(IJ))
	Y1(IJ)=(PRAUT(IJ)+PRACW(IJ))*D2T

c       QC(IJ)=QC(IJ)-Y1(IJ)
c        IF (QC(IJ).LT.0.0) THEN
c            IF(Y1(IJ) .EQ. 0.) Y1(IJ)=CMAX2
c          A1=QC(IJ)/Y1(IJ)+1.
c         PRAUT(IJ)=PRAUT(IJ)*A1
c         PRACW(IJ)=PRACW(IJ)*A1
c         QC(IJ)=0.0
c        ENDIF

	if (qc(ij) .lt. y1(ij) .and. y1(ij) .ge. cmin2) then
	    y2(ij)=qc(ij)/(y1(ij)+cmin2)
	   praut(ij)=praut(ij)*y2(ij)
	   pracw(ij)=pracw(ij)*y2(ij)
	   qc(ij)=0.0
	 else
	   qc(ij)=qc(ij)-y1(ij)
	 endif
c
	 PR(IJ)=(PRAUT(IJ)+PRACW(IJ))*D2T
	 QR(IJ)=QR(IJ)+PR(IJ)
C     ******   CONDENSATION   ******
      CND(IJ)=DM(IJ)/(1.+RSUB1(IJ))
      DD(IJ)=-CND(IJ)-QC(IJ)
      IF(DD(IJ) .LE. 0.0) DD(IJ)=0.0
	DD1(IJ)=-QC(IJ)
      IF(CND(IJ) .LE. DD1(IJ)) CND(IJ)=DD1(IJ)
C     ******   EVAPORATION   ******
      RQ(IJ)=R00*QR(IJ)
CCCCC
      ERN(IJ)=(((1.6+124.9*RQ(IJ)**.2046)*RQ(IJ)**.525)
     1     /(2.55E6/(P00*QVS(IJ))+5.4E5))*(-DM(IJ)/(R00*QVS(IJ)))*D2T
      IF(ERN(IJ) .LE. 0.0) ERN(IJ)=0.0
      IF(ERN(IJ) .GE. DD(IJ)) ERN(IJ)=DD(IJ)
      IF(ERN(IJ) .GE. QR(IJ)) ERN(IJ)=QR(IJ)
CC     ******   SATURATION ADJUSTMENT   ******
      QC(IJ)=QC(IJ)+CND(IJ)
      QR(IJ)=QR(IJ)-ERN(IJ)
      PT(IJ)=PT(IJ)+F00*(CND(IJ)-ERN(IJ))
      QV(IJ)=QV(IJ)-CND(IJ)+ERN(IJ)
C      HUM(IJ,K)=(DQV(IJ,K)+Q00)/(QVS(IJ)-RSUB1(IJ)*(ERN(IJ)-CND(IJ)))
CC
CCCCC
CC
C        IF (QV(IJ)+QB0 .LE. 0.) QV(IJ)=-QB0
	IF (QC(IJ) .LE. CMIN) QC(IJ)=0.
	IF (QR(IJ) .LE. CMIN) QR(IJ)=0.
C        Y1(IJ)=PT(IJ)-XX0(IJ)
C        Y2(IJ)=QV(IJ)-XX00(IJ)
C       DPT(IJ,K)=PT(IJ)+TMD(IJ,K)
C       DQV(IJ,K)=QV(IJ)+QMD(IJ,K)
C        TMD(IJ,K)=Y1(IJ)
C        QMD(IJ,K)=Y2(IJ)
	DPT(IJ,K)=PT(IJ)
	DQV(IJ,K)=QV(IJ)
	QCL(IJ,K)=QC(IJ)
	QRN(IJ,K)=QR(IJ)
  250   CONTINUE
CC*********************************************************************
       SCC=0.
       SEE=0.
       DO IJ=ISTART,IJLES
	 DD(IJ)=MAX(-CND(IJ), 0.)
	 CND(IJ)=MAX(CND(IJ), 0.)
	 SCC=SCC+CND(IJ)
	 SEE=SEE+DD(IJ)+ERN(IJ)
	ENDDO
	SC(K)=SCC+SC(K)
	SE(K)=SEE+SE(K)
CC    ***   STATISTICS FOR CONVECTIVE AND ANVIL REGIMES   ***********
       IF (ID .EQ. 1) THEN
	 RDTS=1./D2T
	DO 280 IJ=ISTART,IJLES
	   CND(IJ)=CND(IJ)*RDTS
	   DD(IJ)=DD(IJ)*RDTS
	   ERN(IJ)=ERN(IJ)*RDTS
C         DDA(IJ)=RSW(IJ,K)
C         DDB(IJ)=RLW(IJ,K)
	 Y1(IJ)=QC(IJ)+QR(IJ)
	 DM(IJ)=A0*(RHO1(K)*WW1(IJ,K)+RHO1(K+1)*WW1(IJ,K+1))
	 RQ(IJ)=.005*(WW1(IJ,K)+WW1(IJ,K+1))
  280    CONTINUE
C       DO 1050 KC=1,7
	KC=4
	 DO 320 MT=1,4
	 DO IJ=ISTART,IJLES
	  IBZ(IJ,MT)=0
	   IF(ICS5(IJ,MT).EQ.1) IBZ(IJ,MT)=1
	  ENDDO
  320   CONTINUE
	DO 315 IJ=ISTART,IJLES
	    IBZ(IJ,1)=1
  315    CONTINUE
	 DO MT=1,4
	 DO 330 IJ=ISTART,IJLES
	   IF(KC.EQ.4) GO TO 330
	   IF(KC.LE.3) GO TO 36
	    IF (RQ(IJ).GT.RBY(KC)) IBZ(IJ,MT)=0
	    GO TO 330
   36       IF (RQ(IJ).LT.RBY(KC)) IBZ(IJ,MT)=0
  330    CONTINUE
	 ENDDO
	 DO 350 MT=1,4
	  SWW=0.0
	  SCC=0.0
	  SEE=0.0
	  A1=0.0
	  A9=0.0
	  A10=0.0
	  A11=0.0
	  A33=0.0
	  A34=0.0
	  A38=0.0
	  A43=0.0
	  A44=0.0
	  A45=0
	  A46=0.
	  A48=0.
	  A49=0.
	  A50=0.
	 DO 30 IJ=ISTART,IJLES
	    ASSS(IJ)=CSTTT(IJ)
	  IF(IBZ(IJ,MT).EQ.1) THEN
	   SWW=SWW+DM(IJ)*ASSS(IJ)
	   SCC=SCC+CND(IJ)*ASSS(IJ)
	   SEE=SEE+DD(IJ)*ASSS(IJ)
	   A9=A9+PRAUT(IJ)*ASSS(IJ)
	   A10=A10+PRACW(IJ)*ASSS(IJ)
	   A33=A33+QC(IJ)*ASSS(IJ)
	   A34=A34+QR(IJ)*ASSS(IJ)
	   A38=A38+ERN(IJ)*ASSS(IJ)
	    A43=A43+(QV(IJ)+QA1(K)-QB(K))*ASSS(IJ)
	    A44=A44+(PT(IJ)+TA1(K)-TB(K))*ASSS(IJ)
	    A45=A45+ASSS(IJ)
	    A46=A46+Y1(IJ)*ASSS(IJ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   A48=A48+(Y1(IJ)-QCL1(IJ,K)-QRN1(IJ,K))*ASSS(IJ)
	   A49=A49+(QV(IJ)-DQV1(IJ,K))*ASSS(IJ)
	   A50=A50+(PT(IJ)-DPT1(IJ,K))*ASSS(IJ)
	  ENDIF
   30   CONTINUE
	SMF0(K,MT,KC)=SWW+SMF0(K,MT,KC)
	COC(K,MT,KC)=SCC+COC(K,MT,KC)
	COE(K,MT,KC)=SEE+COE(K,MT,KC)
	RAUT(K,MT,KC)=RAUT(K,MT,KC)+A9
	RACW(K,MT,KC)=RACW(K,MT,KC)+A10
	ERNS(K,MT,KC)=ERNS(K,MT,KC)+A38
	QC0(K,MT,KC)=A33
	QR0(K,MT,KC)=A34
	QV0(K,MT,KC)=A43
	TT0(K,MT,KC)=A44
	SGPT(K,MT,KC)=A45
	TSQQ(K,MT,KC)=A46
	SQHDT(K,MT,KC)=SQHDT(K,MT,KC)+A48
	SQVDT(K,MT,KC)=SQVDT(K,MT,KC)+A49
	SQTDT(K,MT,KC)=SQTDT(K,MT,KC)+A50
	SQC0(K,MT,KC)=SQC0(K,MT,KC)+QC0(K,MT,KC)
	SQR0(K,MT,KC)=SQR0(K,MT,KC)+QR0(K,MT,KC)
	SQV0(K,MT,KC)=SQV0(K,MT,KC)+QV0(K,MT,KC)
	STT0(K,MT,KC)=STT0(K,MT,KC)+TT0(K,MT,KC)
	SSGPT(K,MT,KC)=SSGPT(K,MT,KC)+SGPT(K,MT,KC)
	TSQQ1(K,MT,KC)=TSQQ1(K,MT,KC)+TSQQ(K,MT,KC)
  350   CONTINUE
       ENDIF
       IF(ID.EQ.1) THEN
       DO 42 IJ=ISTART,IJLES
	CND(IJ)=CND(IJ)*RFT*ASSS(IJ)
	ERN(IJ)=(ERN(IJ)+DD(IJ))*RFT*ASSS(IJ)
	Y5(IJ)=DM(IJ)*ASSS(IJ)
	A1=1.E6*R00*QR(IJ)                                                  
	A11=UCOR*(MAX(1.E-5,A1))**1.75
	ZDRY=MAX(1.,A11)
	DBZ(IJ)=10.*LOG10(ZDRY)
	DBZ(IJ)=DBZ(IJ)*ASSS(IJ)
	QC(IJ)=QC(IJ)*ASSS(IJ)
	QR(IJ)=QR(IJ)*ASSS(IJ)
	 iww=0
	 jww=0
	 if(rq(ij) .gt. -0.5) then
	  iww=min(rq(ij)+0.5, 15.) + 1
	 else
	  jww=max(rq(ij)-0.5, -5.0)
	  jww=iabs(jww)
	 endif
   42  CONTINUE
       DO 44 IJ=ISTART,IJLES
	IF(RQ(IJ) .GE. 0.) SCU1(K)=SCU1(K)+CND(IJ)
	IF(RQ(IJ) .LT. 0.) SED1(K)=SED1(K)+ERN(IJ)
   44  CONTINUE
       DO 40 IJ=ISTART,IJLES
	 iww=0
	 jww=0
	 if(rq(ij) .gt. -0.5) then
	  iww=min(rq(ij)+0.5, 15.) + 1
	 else
	  jww=max(rq(ij)-0.5, -5.0)
	  jww=iabs(jww)
	 endif
	IF (ICS(IJ,2) .EQ. 1 .AND. IWW.GE.1) THEN
	  S9(IWW,K)=S9(IWW,K)+CND(IJ)
	  S10(IWW,K)=S10(IWW,K)+ERN(IJ)
	  S15(IWW,K)=S15(IWW,K)+Y5(IJ)
	  S16(IWW,K)=S16(IWW,K)+QC(IJ)
	  S17(IWW,K)=S17(IWW,K)+QR(IJ)
	  S21(IWW,K)=S21(IWW,K)+DBZ(IJ)
	  SCNT(IWW,K)=SCNT(IWW,K)+ASSS(IJ)
	ENDIF
	IF (ICS(IJ,2) .EQ. 1 .AND. jww.GE.1) THEN
	  SN9(jww,K)=SN9(jww,K)+CND(IJ)
	  SN10(jww,K)=SN10(jww,K)+ERN(IJ)
	  SN15(jww,K)=SN15(jww,K)+Y5(IJ)
	  SN16(jww,K)=SN16(jww,K)+QC(IJ)
	  SN17(jww,K)=SN17(jww,K)+QR(IJ)
	  SN21(jww,K)=SN21(jww,K)+DBZ(IJ)
	  SNCNT(jww,K)=SNCNT(jww,K)+ASSS(IJ)
	ENDIF
	IF (ICS(IJ,3) .EQ. 1 .AND. IWW.GE.1) THEN
	  SS9(IWW,K)=SS9(IWW,K)+CND(IJ)
	  SS10(IWW,K)=SS10(IWW,K)+ERN(IJ)
	  SS15(IWW,K)=SS15(IWW,K)+Y5(IJ)
	  SS16(IWW,K)=SS16(IWW,K)+QC(IJ)
	  SS17(IWW,K)=SS17(IWW,K)+QR(IJ)
	  SS21(IWW,K)=SS21(IWW,K)+DBZ(IJ)
	  SSCNT(IWW,K)=SSCNT(IWW,K)+ASSS(IJ)
	ENDIF
	IF (ICS(IJ,3) .EQ. 1 .AND. jww.GE.1) THEN
	  SSN9(jww,K)=SSN9(jww,K)+CND(IJ)
	  SSN10(jww,K)=SSN10(jww,K)+ERN(IJ)
	  SSN15(jww,K)=SSN15(jww,K)+Y5(IJ)
	  SSN16(jww,K)=SSN16(jww,K)+QC(IJ)
	  SSN17(jww,K)=SSN17(jww,K)+QR(IJ)
	  SSN21(jww,K)=SSN21(jww,K)+DBZ(IJ)
	  SSNCNT(jww,K)=SSNCNT(jww,K)+ASSS(IJ)
	ENDIF
   40  CONTINUE
       ENDIF
 1000 CONTINUE
C     ****************************************************************
      IF (ID .EQ. 0) THEN
       LC=0
       LS=0
       DO 390 IJ=ISTART,IJLES
	 LC=LC+ICS(IJ,2)
	 LS=LS+ICS(IJ,3)
  390  CONTINUE
	 IF (LC .EQ. 0) LC=1000000
	 IF (LS .EQ. 0) LS=1000000
       DO 400 MT=1,3
	A1=RIJL2
	IF (MT .EQ. 2) A1=1./FLOAT(LC)
	IF (MT .EQ. 3) A1=1./FLOAT(LS)
       DO K=KT1,KT2
	 B1(K)=0.0
	 B2(K)=0.0
       ENDDO
       DO IJ=ISTART,IJLES
       DO 415 K=KT1,KT2
	IF(ICS(IJ,MT) .EQ. 0) GO TO 415
	 B1(K)=B1(K)+DPT(IJ,K)
	 B2(K)=B2(K)+DQV(IJ,K)
  415  CONTINUE
       ENDDO
       DO 430 K=KT1,KT2
	TB00(K,MT)=B1(K)*A1
	QB00(K,MT)=B2(K)*A1
  430  CONTINUE
  400  CONTINUE
       ENDIF
C     ****************************************************************
      D2T=D22T
      CALL SATDT
      RETURN
      END

      SUBROUTINE SATDT
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880)
      PARAMETER (NZ2=2*NZ,NZ3=3*NZ,NZ4=4*NZ,NM2=2*NM,NM11=11*NM)
      PARAMETER (NB3=NZ+NX,NB4=5*NT)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC(4)
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(9),EPS,PSFC(5)

      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B2S/ QCS1(NX,NY,NZ)

      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),ZX(NZ4),AM(NZ),ZQ(NZ3),WB(NZ),ZW(NZ2),RRHO(NZ),WBX(NB3)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      COMMON/BRH1/ SRRO(NZ),QRRO(NZ),SQC(NZ),SQR(NZ),SQI(NZ),SQS(NZ),
     1   SQG(NZ),STQC(NZ),STQR(NZ),STQI(NZ),STQS(NZ),STQG(NZ),TTTD(NB4)
      COMMON/BBA/ Y1(NX,NY),Y2(NX,NY),Y3(NX,NY),Y4(NX,NY),Y7(NM2),
     1  B0(NM),B1(NM),B2(NM),Y6(NM11)
      save
C     ****************************************************************

c$doacross local(k,j,i)
      DO K=2,KLES
	SQC(K)=0.0
	SQR(K)=0.0
	SQI(K)=0.0
	SQS(K)=0.0
	SQG(K)=0.0
	B1(K)=0.0
	B2(K)=0.0
	SQ(K)=0.0
      DO J=2,JLES
      DO I=2,ILES
       SQC(K)=SQC(K)+QCL(I,J,K)
       SQR(K)=SQR(K)+QRN(I,J,K)
       SQI(K)=SQI(K)+QCI(I,J,K)
       SQS(K)=SQS(K)+QCS(I,J,K)
       SQG(K)=SQG(K)+QCG(I,J,K)
      B1(K)=B1(K)+DPT(I,J,K)
      B2(K)=B2(K)+DQV(I,J,K)
      SQ(K)=SQ(K)+QCL(I,J,K)+QRN(I,J,K)+QCI(I,J,K)+QCS(I,J,K)+QCG(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      DO 460 K=2,KLES
       STQC(K)=STQC(K)+SQC(K)*RIJL2
       STQR(K)=STQR(K)+SQR(K)*RIJL2
       STQI(K)=STQI(K)+SQI(K)*RIJL2
       STQS(K)=STQS(K)+SQS(K)*RIJL2
       STQG(K)=STQG(K)+SQG(K)*RIJL2
      B1(K)=B1(K)*RIJL2
      B2(K)=B2(K)*RIJL2
      SQ(K)=SQ(K)*RIJL2
  460 CONTINUE

C     **************************************
      DO 500 K=2,KLES
       B1K=B1(K)
       B2K=B2(K)
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
	DPT(I,J,K)=DPT(I,J,K)-B1K
	DQV(I,J,K)=DQV(I,J,K)-B2K
       ENDDO
       ENDDO
       IF (N.GT.1) GO TO 515
	 IF (IJKADV .EQ. 1) THEN
	   Y6(K)=TA(K)
	   Y7(K)=QA(K)
	  TA(K)=TA(K)+b1(K)
	  QA(K)=QA(K)+b2(K)
          TA1(K)=ta(K)
          QA1(K)=qa(K)
c	   TA1(K)=y6(K)
c	   QA1(K)=y7(K)
	 ELSE
	  TA(K)=TA1(K)+b1(K)
	  QA(K)=QA1(K)+b2(K)
	 ENDIF
       GO TO 516
  515  CONTINUE

	if (ijkadv .eq. 0) then
c$doacross local(j,i)
	  do j=2,jles
	  do i=2,iles
	    dpt1(i,j,k)=dpt1(i,j,k)+eps*dpt(i,j,k)
	    dqv1(i,j,k)=dqv1(i,j,k)+eps*dqv(i,j,k)
	    qcl1(i,j,k)=qcl1(i,j,k)+eps*qcl(i,j,k)
	    qrn1(i,j,k)=qrn1(i,j,k)+eps*qrn(i,j,k)
	  enddo
	  enddo
	endif

	if (iwater .eq. 0 .and. ijkadv .eq. 0) then
c$doacross local(j,i)
	   do j=2,jles
	   do i=2,iles
	     qci1(i,j,k)=qci1(i,j,k)+eps*qci(i,j,k)
	     qcs1(i,j,k)=qcs1(i,j,k)+eps*qcs(i,j,k)
	     qcg1(i,j,k)=qcg1(i,j,k)+eps*qcg(i,j,k)
	   enddo
	   enddo
	endif

        if (ijkadv .eq. 1) then
c$doacross local(j,i)
          do j=2,jles
          do i=2,iles
            dpt1(i,j,k)=dpt(i,j,k)
            dqv1(i,j,k)=dqv(i,j,k)
            qcl1(i,j,k)=qcl(i,j,k)
            qrn1(i,j,k)=qrn(i,j,k)
          enddo
          enddo
        endif

        if (iwater .eq. 0 .and. ijkadv .eq. 1) then
c$doacross local(j,i)
           do j=2,jles
           do i=2,iles
             qci1(i,j,k)=qci(i,j,k)
             qcs1(i,j,k)=qcs(i,j,k)
             qcg1(i,j,k)=qcg(i,j,k)
           enddo
           enddo
        endif

	 IF (IJKADV .EQ. 1) THEN
	   Y6(K)=TA(K)
	   Y7(K)=QA(K)
	  TA(K)=TA(K)+b1(K)
	  QA(K)=QA(K)+b2(K)
           TA1(K)=ta(K)
           QA1(K)=qa(K)
c	    TA1(K)=y6(K)
c	    QA1(K)=y7(K)
	 ELSE
	   A0=TA(K)+EPS*(-2.*TA(K)+TA1(K))
	  TA(K)=TA1(K)+b1(K)
	  TA1(K)=A0+EPS*TA(K)
	   A0=QA(K)+EPS*(-2.*QA(K)+QA1(K))
	  QA(K)=QA1(K)+b2(K)
	  QA1(K)=A0+EPS*QA(K)
	 ENDIF
  516  CONTINUE
       ST(K)=TA(K)-TB(K)
       SV(K)=QA(K)-QB(K)
  500 CONTINUE
C   *************************************************
      CALL BNDOP (DPT,DPT1)
      CALL BNDOP (DQV,DQV1)
      CALL BNDOP (QCL,QCL1)
      CALL BNDOP (QRN,QRN1)
      IF (IWATER .EQ. 0) THEN
	CALL BNDOP (QCI,QCI1)
	CALL BNDOP (QCS,QCS1)
	CALL BNDOP (QCG,QCG1)
      ENDIF
      SQ(1)=SQ(2)
      SQ(KMAX)=SQ(KLES)
       TA(1)=TA(2)
       QA(1)=QA(2)
       TA1(1)=TA1(2)
       QA1(1)=QA1(2)
       TA(KMAX)=2.*TA(KLES)-TA(KL2)
       QA(KMAX)=2.*QA(KLES)-QA(KL2)
       TA1(KMAX)=2.*TA1(KLES)-TA1(KL2)
       QA1(KMAX)=2.*QA1(KLES)-QA1(KL2)
      DO 600 K=3,KLES
	if (ijkadv .eq. 1) then
	  FD(K)=(TA(K)-TA(K-1))*RHO1(K)*.5
	  FE(K)=(QA(K)-QA(K-1))*RHO1(K)*.5
	else
	  FD(K)=(TA(K)-TA(K-1))*RHO1(K)*.5
	  FE(K)=(QA(K)-QA(K-1))*RHO1(K)*.5
	endif
  600 continue 
cc
       FD(KMAX)=FD(KLES)
       FE(KMAX)=FE(KLES)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DOMAIN (IFLAG)
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      PARAMETER (NB=NX*NY+16*NM,NB1=13*NZ+NX)
      PARAMETER (NZ4=4*NZ,NZ2=2*NZ,NZ7=7*NZ)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NIS(3),
     1  K1,K2
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/BSAT/ X(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B5/ TB(NZ4),TA(NZ),QA(NZ),TA1(NZ),QA1(NZ),COEF(NB1)
      COMMON/B6/ FD(NZ2),P0(NZ),PI(NZ),F0(NZ7)
      COMMON/BBA/ TAIR(NX,NY),QSS(NX,NY),QV(NX,NY),XYPM(NB)
      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	K1=2
	K2=KLES
       IF(IFLAG .EQ. 1) RETURN
       T0=273.16

      DO K=2,KLES
	RP0=C38/P0(K)
	TAK=TA(K)
	QAK=QA(K)
	PIK=PI(K)
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
	TAIR(I,J)=(DPT(I,J,K)+TAK)*PIK
	IF(TAIR(I,J) .GE. T0) THEN
	 QSS(I,J)=RP0*EXP(C172-C409/(TAIR(I,J)-C358))
	ELSE
	 QSS(I,J)=RP0*EXP(C218-C580/(TAIR(I,J)-C76))
	ENDIF
	QV(I,J)=MAX(0., DQV(I,J,K)+QAK-QSS(I,J))
       X(I,J,K)=QV(I,J)+QCL(I,J,K)+QCI(I,J,K)+QRN(I,J,K)
     1                 +QCS(I,J,K)+QCG(I,J,K)
       ENDDO
       ENDDO
       ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Y=0.
      DO 300 K=2,KLES
c$doacross local(j,i), reduction(y)
       DO J=2,JLES
       DO I=2,ILES
	  Y=MAX(Y,X(I,J,K))
       ENDDO
       ENDDO
       IF(Y .GT. 0.) THEN
	K1=K
	GO TO 35
       ENDIF
  300 CONTINUE
   35 Y=0.
      DO 400 K=KLES,2,-1
c$doacross local(j,i), reduction(y)
       DO J=2,JLES
       DO I=2,ILES
	  Y=MAX(Y,X(I,J,K))
       ENDDO
       ENDDO
       IF(Y .GT. 0.) THEN
	K2=K
	GO TO 45
       ENDIF
  400 CONTINUE
   45 K1=MAX(2,K1-2)
      K2=MIN(KLES,K2+2)
      RETURN
      END

      SUBROUTINE DOMAINW (IFLAG)
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      PARAMETER (NB=NX*NY+16*NM,NB1=13*NZ+NX)
      PARAMETER (NZ4=4*NZ,NZ2=2*NZ,NZ7=7*NZ)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NIS(3),
     1  K1,K2
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/BSAT/ X(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B5/ TB(NZ4),TA(NZ),QA(NZ),TA1(NZ),QA1(NZ),COEF(NB1)
      COMMON/B6/ FD(NZ2),P0(NZ),PI(NZ),F0(NZ7)
      COMMON/BBA/ TAIR(NX,NY),QSS(NX,NY),QV(NX,NY),XYPM(NB)
      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	K1=2
	K2=KLES
       IF(IFLAG .EQ. 1) RETURN
      DO K=2,KLES
	RP0=C38/P0(K)
	TAK=TA(K)
	QAK=QA(K)
	PIK=PI(K)
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
	TAIR(I,J)=(DPT(I,J,K)+TAK)*PIK
	 QSS(I,J)=RP0*EXP(C172-C409/(TAIR(I,J)-C358))
	QV(I,J)=MAX(0., DQV(I,J,K)+QAK-QSS(I,J))
       X(I,J,K)=QV(I,J)+QCL(I,J,K)+QRN(I,J,K)
       ENDDO
       ENDDO
       ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Y=0.
      DO 300 K=2,KLES
c$doacross local(j,i), reduction(y)
       DO J=2,JLES
       DO I=2,ILES
	  Y=MAX(Y,X(I,J,K))
       ENDDO
       ENDDO
       IF(Y .GT. 0.) THEN
	K1=K
	GO TO 35
       ENDIF
  300 CONTINUE
   35 Y=0.
      DO 400 K=KLES,2,-1
c$doacross local(j,i), reduction(y)
       DO J=2,JLES
       DO I=2,ILES
	  Y=MAX(Y,X(I,J,K))
       ENDDO
       ENDDO
       IF(Y .GT. 0.) THEN
	K2=K
	GO TO 45
       ENDIF
  400 CONTINUE
   45 K1=MAX(2,K1-2)
      K2=MIN(KLES,K2+2)
      RETURN
      END

      SUBROUTINE AKCOEF
C     ****   COMPUTE DIFFUSION COEFFICIENT
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      PARAMETER (NM12=12*NM)
c      COMMON/RTIME/ ISEC1
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1  NRAN,KT1,KT2
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(4),AL,CP,RA,CK,CE,EPS,PSFC(5)
C
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSKT2(NZ),BSK(NZ),BSK4(NZ),
     1   BSIT(NX,NY),BSIT2(NX,NY),BSI(NX,NY),BSI4(NX,NY),
     2   BSJT(NX,NY),BSJT2(NX,NY),BSJ(NX,NY),BSJ4(NX,NY)
      COMMON/B4Z/ BSITZ(NZ),BSIT2Z(NZ),BSIZ(NZ),BSI4Z(NZ),BSJTZ(NZ),
     1   BSJT2Z(NZ),BSJZ(NZ),BSJ4Z(NZ)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
C
      COMMON/BADV/ U1(NX,NY,NZ)
      COMMON/BADV1/ UU1(NX,NY,NZ)
      COMMON/BSAT/ VV1(NX,NY,NZ)
      COMMON/BSAT1/ WW1(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
C      COMMON/B1TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
C      COMMON/B1CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
C      COMMON/B1IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B4UV/ UUD(NX,NY,NZ),VVD(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
      COMMON/B4WP/ WWD(NX,NY,NZ)
      COMMON/B4A/ AK1(NX,NY,NZ)
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
C
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
C     COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA1(NZ),QA1(NZ),TA(NZ),
C    1  QA(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
C    2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
C
      COMMON/BBA/ TP(NX,NY),TM(NX,NY),QP(NX,NY),QM(NX,NY),Y1(NM),Y2(NM),
     1   Y3(NM),Y4(NM),Y12(NM12)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
      DIMENSION XY1(NX,NY),XY2(NX,NY),XY3(NX,NY),XY4(NX,NY),XY5(NX,NY),
     1   XY6(NX,NY),XY7(NX,NY)

      save

C   ********************************************************************
       a0000=1. !!sh
       if (imlifting .eq. 0) a0000=0. !!sh

       a_irf=1.
	 if (irf .eq. 0) a_irf=0.
c
       awks=1.
	 if (ijkadv .eq. 0) awks=1.
c
      ALCP=AL/CP
      CPRA=CP*RA
      ALSQ=.622*AL*AL
      CPAL=.622*1.61*CP*AL
	C13=1./2.
      DO K=1,KMAX
       Y1(K)=1./TA1(K)
       Y2(K)=1./PI(K)
       Y3(K)=-AM(K)*Z1(K)*RD2Z
      ENDDO
      if (ijkadv .eq. 1) then
c$doacross local(k,j,i)
	 DO K=1,KMAX
	 DO J=1,JMAX
	 DO I=1,IMAX
	    UU1(I,J,K)=.5*(3.*UUD(I,J,K)-UMD(I,J,K))
	    VV1(I,J,K)=.5*(3.*VVD(I,J,K)-VMD(I,J,K))
	    WW1(I,J,K)=.5*(3.*WWD(I,J,K)-WMD(I,J,K))
	 ENDDO
	 ENDDO
	 ENDDO
      else
c$doacross local(k,j,i)
	 DO K=1,KMAX
	 DO J=1,JMAX
	 DO I=1,IMAX
	    UU1(I,J,K)=UUD(I,J,K)
	    VV1(I,J,K)=VVD(I,J,K)
	    WW1(I,J,K)=WWD(I,J,K)
	 ENDDO
	 ENDDO
	 ENDDO
      endif
c$doacross local(k,j,i)
      DO K=1,KMAX
      DO J=1,JMAX
      DO I=1,IMAX
	 AK(I,J,K)=C13*AK(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C    ****************************************
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
       TM(I,J)=TA1(1)+DPT1(I,J,1)
       QM(I,J)=QA1(1)+DQV1(I,J,1)
       XY1(I,J)=TA1(2)+DPT1(I,J,2)
	XY2(I,J)=QA1(2)+DQV1(I,J,2)
      ENDDO
      ENDDO
      DO 2001 K=2,KLES
	KP=K+1
	KM=K-1
	TA1KP=TA1(KP)
	QA1KP=QA1(KP)
	Y1K=Y1(K)
	PIK=PI(K)
	Y2KP=Y2(KP)
	Y2KM=Y2(KM)
	Y3K=Y3(K)
	TA1K=TA1(K)

c$doacross local(j,i)
      DO 2000 J=2,JLES
       DO I=2,ILES
	 XY4(I,J)=QCL1(I,J,K)+QCI1(I,J,K)
	 TP(I,J)=TA1KP+DPT1(I,J,KP)
	 QP(I,J)=QA1KP+DQV1(I,J,KP)
	IF (XY4(I,J) .LT. 1.E-5) THEN
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 U1(I,J,K)=Y1K*(TP(I,J)*(1.+.61*QP(I,J))
     1                   -TM(I,J)*(1.+.61*QM(I,J)))
	ELSE
	  XY3(I,J)=XY1(I,J)*PIK
	  XY4(I,J)=CPRA*XY3(I,J)*XY3(I,J)
	  XY5(I,J)=(XY4(I,J)+CPAL*XY3(I,J)*XY2(I,J))
     1                /((XY4(I,J)+ALSQ*XY2(I,J))*TA1K)
	 XY6(I,J)=TP(I,J)-TM(I,J)+ALCP*(QP(I,J)*Y2KP-QM(I,J)*Y2KM)
	  XY7(I,J)=QP(I,J)-QM(I,J)+QCL1(I,J,KP)
     1                  -QCL1(I,J,KM)+QCI1(I,J,KP)-QCI1(I,J,KM)
	 U1(I,J,K)=XY5(I,J)*XY6(I,J)-XY7(I,J)
	ENDIF
       U1(I,J,K)=Y3K*U1(I,J,K)
c       ENDDO
c        DO I=2,ILES
	 TM(I,J)=XY1(I,J)
	 QM(I,J)=XY2(I,J)
	 XY1(I,J)=TP(I,J)
	 XY2(I,J)=QP(I,J)
       ENDDO
 2000  CONTINUE

 2001  CONTINUE
c
c$doacross local(j,i)
      DO J=1,JMAX
      DO I=1,IMAX
       QM(I,J)=AK1(I,J,1)*AK1(I,J,1)
       TM(I,J)=AK1(I,J,2)*AK1(I,J,2)
      ENDDO
      ENDDO


      DO 3000 K=2,KLES
	KP=K+1
	KM=K-1
	AMK=AM(K)
	AM1K=AM1(K)
	AM1KP=AM1(KP)
	a22z=AMK*RRHO(K)*RD2Z
c        amm2=amk*amk*R4DZ2
	ammkp=awks*am1kp*rho1(kp)*bsk(kp)
	ammk=awks*am1k*rho1(k)*bsk(k)
	A33K=RD4Z*AMK
	A44K=AMK*AMK*RDZ2
	A55K=AMK*R2DZ2
       COEFK=COEF(K)
       WBK=WB(K)
       WBKP=WB(KP)
       Z2K=Z2(K)

c$doacross local(j,i)
       DO J=1,JMAX
       DO I=1,IMAX
	 TP(I,J)=AK1(I,J,KP)*AK1(I,J,KP)
       ENDDO
       ENDDO

c$doacross local(j,jp,jm,i,ip,im)
       DO J=2,JLES
       JP=J+1
       JM=J-1
       DO I=2,ILES
       IP=I+1
       IM=I-1
       XY1(I,J)=
     1  (WW1(IP,J,KP)+WW1(IP,J,K)-WW1(IM,J,KP)-WW1(IM,J,K))*RD4X
     2  +(UU1(IP,J,KP)+UU1(I,J,KP)-UU1(IP,J,KM)-UU1(I,J,KM))*A33K
	XY1(I,J)=XY1(I,J)*XY1(I,J)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       XY2(I,J)=
     1  (UU1(IP,JP,K)+UU1(I,JP,K)-UU1(IP,JM,K)-UU1(I,JM,K))*RD4Y
     2  +(VV1(IP,J,K)+VV1(IP,JP,K)-VV1(IM,J,K)-VV1(IM,JP,K))*RD4X
	XY2(I,J)=XY2(I,J)*XY2(I,J)
       XY3(I,J)=
     1  (WW1(I,JP,KP)+WW1(I,JP,K)-WW1(I,JM,K)-WW1(I,JM,KP))*RD4Y
     2  +(VV1(I,J,KP)+VV1(I,JP,KP)-VV1(I,J,KM)-VV1(I,JP,KM))
     3  *A33K
	XY3(I,J)=XY3(I,J)*XY3(I,J)
       QP(I,J)=COEFK*(XY1(I,J)+XY2(I,J)+XY3(I,J)
     1  +((UU1(IP,J,K)-UU1(I,J,K))**2*RDX2
     2  +(VV1(I,JP,K)-VV1(I,J,K))**2*RDY2
     3  +A44K*(WW1(I,J,KP)-WW1(I,J,K)
     4         +a0000*Y0(I,J)*(WBKP-WBK))**2 )*2.)
       U1(I,J,K)=U1(I,J,K)+QP(I,J)-Z2K*TM(I,J)
     1  +((TM(IP,J)-TM(I,J))-(TM(I,J)-TM(IM,J)))*R2DX2
     2  +((TM(I,JP)-TM(I,J))-(TM(I,J)-TM(I,JM)))*R2DY2
     3  +A55K*(AM1KP*(TP(I,J)-TM(I,J))-AM1K*(TM(I,J)-QM(I,J)))
       ENDDO
       ENDDO

c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
	U1(I,J,K)=U1(I,J,K)-A22Z*(
     1      -2.*(AMMKP*(AK1(I,J,KP)-AK1(I,J,K))
     2      -AMMK*(AK1(I,J,K)-AK1(I,J,KM)))*RDZ)-a0000*Y0(I,J)
     1      *AM(K)*(WB(KP)+WB(K))*(AK(I,J,KP)-AK(I,J,KM))*RD4Z
       ENDDO
       ENDDO

c$doacross local(j,i,ipp)
       DO J=2,JLES
	 do i=2,iles
	   ipp=i+2
	     if (i.eq.iles) ipp=3
	xy2(i,j)=a4k(k)*(ak1(ipp,j,k)-3.*(ak1(i+1,j,k)-ak1(i,j,k))
     1                   -ak1(i-1,j,k))
	enddo
	enddo

c$doacross local(j,i,im)
       DO J=2,JLES
	 do i=2,iles
	 im=i-1
	 if (i.eq.2) im=iles
	 u1(i,j,k)=u1(i,j,k)+awks*(xy2(im,j)-xy2(i,j))
	enddo
	enddo

c$doacross local(j,jp,jm,jpp,i)
       do j=2,jles
	 jp=j+1
	 jm=j-1
	 jpp=j+2
	 if(j.eq.jles) jpp=3
	 DO I=2,ILES
	 xy1(i,j)=a4k(k)*(ak1(i,jpp,k)-3.*(ak1(i,jp,k)-ak1(i,j,k))
     1                   -ak1(i,jm,k))
	enddo
	enddo

c$doacross local(j,jm,i)
       do j=2,jles
	 jm=j-1
	 if (j.eq.2) jm=jles
	 DO I=2,ILES
	    u1(i,j,k)=u1(i,j,k)+awks*(xy1(i,jm)-xy1(i,j))
	enddo
	enddo

c$doacross local(j,i)
	DO J=1,JMAX
	DO I=1,IMAX
	 QM(I,J)=TM(I,J)
	TM(I,J)=TP(I,J)
	ENDDO
        ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        i_four=0
c
        if (i_four .eq.10) then
c$doacross local(j,i,im)
	do j=2,jles
	  do i=2,imax
	   im=i-1
	   xy1(i,j)=-(ak(i,j,k)+ak(im,j,k))*
     1            (ak1(i,j,k)-ak1(im,j,k))*r2dx2
	  enddo
	enddo

c$doacross local(j,i)
	do j=2,jles
	  do i=2,iles
	   u1(i,j,k)=u1(i,j,k)+xy1(i,j)-xy1(i+1,j)
	  enddo
	enddo

c$doacross local(j,jm,i)
	do j=2,jmax
	 jm=j-1
	do i=2,iles
	   xy1(i,j)=-(ak(i,j,k)+ak(i,jm,k))*
     1            (ak1(i,j,k)-ak1(i,jm,k))*r2dy2
	enddo
	enddo

c$doacross local(j,i)
	do j=2,jles
	do i=2,iles
	   u1(i,j,k)=u1(i,j,k)+xy1(i,j)-xy1(i,j+1)
	enddo
	enddo
c
	endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 3000  CONTINUE

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (i_four .eq.10) then
c
c$doacross local(j,i)
       do j=2,jles
       do i=2,iles
	  qm(i,j)=0.
       enddo
       enddo

      do 3600 k=3,kmax
	  km=k-1
	  am1k=am1(k)
	  rho1k=rho1(k)
c$doacross local(j,i)
	 do j=2,jles
	 do i=2,iles
c          tm(i,j)=rho1(k)*(-(am1(k)*(ak1(i,j,k)-ak1(i,j,km))
c     1                  *(ak(i,j,k)+ak(i,j,km)))*rdz)
	  tm(i,j)=rho1k*(-(am1k*(ak1(i,j,k)-ak1(i,j,km))
     1                  *(ak(i,j,k)+ak(i,j,km)))*rdz)
	 enddo
	 enddo

	 if (km .eq. kles) then
c$doacross local(j,i)
	   do j=2,jles
	   do i=2,iles
	     tm(i,j)=0.
	   enddo
	   enddo
       endif

c$doacross local(j,i)
       do j=2,jles
       do i=2,iles
	u1(i,j,km)=u1(i,j,km)+am(km)*(qm(i,j)-tm(i,j))*rd2z*rrho(km)
	 qm(i,j)=tm(i,j)
       enddo
       enddo
 3600 continue

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     ******************************************************************
      if (ijkadv .eq. 1) then
C  
c$doacross local(k,rfak,j,i)
	 do k=2,kles
	    rfak=a_irf*rfa(k)
	 do j=2,jles
	 do i=2,iles
	    U1(I,J,K)=U1(I,J,K)-RFAK*AK1(I,J,K)
	 enddo
	 enddo
	 enddo

       CALL ADVECTAK (AK,AK1)
cc
ccc
      else
ccc
cc
	DO 4000 K=2,KLES
	  kp=k+1
	  km=k-1
	    a11k=.5*rrho(k)*am(k)*rdz
	  DO J=2,JLES
	  DO I=2,ILES
	     u1(i,j,k)=u1(i,j,k)-a11k*
     1                  (rho1(kp)*ww1(i,j,kp)*(ak(i,j,k)+ak(i,j,kp))
     2                   -rho1(k)*ww1(i,j,k)*(ak(i,j,k)+ak(i,j,km)))
     3                          -a_irf*RFA(K)*AK1(I,J,K)
c
	  ENDDO
	  ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	  IF(IADVH.EQ.2)THEN
c$doacross local(j,i)
	    do j=2,jles
	    do i=2,imax
	    xy1(i,j)=uu1(i,j,k)*(ak(i,j,k)+ak(i-1,j,k))*rd2x
	    enddo
	    enddo
c$doacross local(j,jm,i)
	    do j=2,jmax
	      jm=j-1
	      do i=2,iles
	    xy3(i,j)=vv1(i,j,k)*(ak(i,j,k)+ak(i,jm,k))*rd2y
	      enddo
	    enddo
c$doacross local(j,jp,i)
	    do j=2,jles
	      jp=j+1
	      do i=2,iles
	      u1(i,j,k)=u1(i,j,k)+(xy1(i,j)-xy1(i+1,j))
     1                         +(xy3(i,j)-xy3(i,jp))
	     enddo
	    enddo
	  ELSE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ******   4-th order  horizontal advection terms   ****************

c$doacross local(j,i,ipp)
	   DO J=2,JLES
	   do i=2,iles
	      ipp=i+2
	      if (i.eq.iles) ipp=3
	      xy2(i,j)=uu1(i+1,j,k)*(d58x*(ak(i,j,k)+ak(i+1,j,k))
     1                -d16x*(ak(i-1,j,k)+ak(ipp,j,k)))
	   enddo
	   enddo

c$doacross local(j,i)
	   DO J=2,JLES
	   do i=2,imax
	     xy3(i,j)=-d48x*uu1(i,j,k)*(ak(i,j,k)+ak(i-1,j,k))
	   enddo
	   enddo

c$doacross local(j,i,im,ipp)
	   DO J=2,JLES
	   do i=2,iles
	      im=i-1
	       if(i.eq.2) im=iles
	      ipp=i+2
	       if (i.eq.iles) ipp=3
	     u1(i,j,k)=u1(i,j,k)+(xy2(im,j)-xy2(i,j)
     1                           +xy3(im,j)-xy3(ipp,j))
	   enddo
	   enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   DO 470 I=2,ILES
	   do j=2,jles
	      jpp=j+2
	       if(j.eq.jles) jpp=3
	     xy1(i,j)=vv1(i,j+1,k)*(d58y*(ak(i,j,k)+ak(i,j+1,k))
     1                            -d16y*(ak(i,j-1,k)+ak(i,jpp,k)))
	   enddo
	   do j=2,jmax
	      xy3(i,j)=-d48y*vv1(i,j,k)*(ak(i,j,k)+ak(i,j-1,k))
	   enddo
	   do j=2,jles
	      jm=j-1
	       if(j.eq.2) jm=jles
	      jpp=j+2
	       if(j.eq.jles) jpp=3
	     u1(i,j,k)=u1(i,j,k)+(xy1(i,jm)-xy1(i,j)
     1                           +xy3(i,jm)-xy3(i,jpp))
	   enddo
  470      CONTINUE
c
c
cc
ccc
	  endif
c
 4000   CONTINUE

	  y1(2)=0.
	  y1(kmax)=0.

	  r2dt8=0.008/d2t

	  do j=2,jles
	  do i=2,iles

	    do k=3,kles
	       y1(k)=ak1(i,j,k)-ak1(i,j,k-1)
	    enddo

	    do k=2,kles
	       u1(i,j,k)=u1(i,j,k)-r2dt8*(y1(k)-y1(k+1))
	    enddo

	  enddo
	  enddo

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
cc
c
cccccc
      DO K=2,KLES
	 A1=TBSK(K)
	if (ijkadv .eq. 0) then
c$doacross local(j,i)
	   DO J=2,JLES
	   DO I=2,ILES
	      XY1(I,J)=AK(I,J,K)+EPS*(-2.*AK(I,J,K)+AK1(I,J,K))
	     AK(I,J,K)=AK1(I,J,K)+U1(I,J,K)*D2T
	     AK1(I,J,K)=XY1(I,j)+EPS*AK(I,J,K)
	  ENDDO
	  ENDDO
	else
c$doacross local(j,i)
	  DO J=2,JLES
	  DO I=2,ILES
	      AK(I,J,K)=AK(I,J,K)+U1(I,J,K)*DT
              AK1(I,J,K)=AK(I,J,K)
	  ENDDO
	  ENDDO
	endif
c$doacross local(j,i)
	DO J=2,JLES
	DO I=2,ILES
	   AK(I,J,K)=MIN(A1, MAX(AK(I,J,K), 0.0))
	   AK1(I,J,K)=MIN(A1, MAX(AK1(I,J,K), 0.0))
	ENDDO
	ENDDO
      ENDDO
C     ***************************************************************
      CALL BNDOP (AK,AK1)
      RETURN
      END

      SUBROUTINE ADVECTAK (X,X1)
CC    ****   COMPUTE ADVECTION OF K COEFFICIENT
      PARAMETER (NX=130,NY=130,NZ=43)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
      DIMENSION X(NX,NY,1),X1(NX,NY,1)
CC    LOCAL VARIABLES
      COMMON /tmp1/ U1(NX,NY,NZ),V1(NX,NY,NZ),W1(NX,NY,NZ)
      save
CC    **************************************************************
c$doacross local(k,j,i)
      DO K=1,KMAX
      DO J=1,JMAX
      DO I=1,IMAX
        U1(I,J,K)=.5*(3.*UU1(I,J,K)-UMD(I,J,K))
        V1(I,J,K)=.5*(3.*VV1(I,J,K)-VMD(I,J,K))      
        W1(I,J,K)=.5*(3.*WW1(I,J,K)-WMD(I,J,K))
      ENDDO
      ENDDO
      ENDDO
c$doacross local(k,j,i)
      DO K=1,KMAX
      DO J=1,JMAX
      DO I=1,IMAX
         IF(X(I,J,K) .LT. 1.E-5) X(I,J,K)=0.
         X1(I,J,K)=X(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C
      CALL FADV (X,U1,V1,W1)
CC    NEED TO COMPUTE BOUNDARY CONDITION
       call bndop (x,x1)
CC    *********************************************
      CALL FADVUW (X,X1)
      CALL FADV (X,U1,V1,W1)
       call bndop (x,x1)
C
c$doacross local(k,j,i)
      DO K=1,KMAX
      DO J=1,JMAX
      DO I=1,IMAX
         IF(X(I,J,K) .LT. 1.E-5) X(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
      SUBROUTINE ADVECTW (SQ)
C     ****   COMPUTE W1 AND FORCE TERMS
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      PARAMETER (NM10=10*NM)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1  NRAN,KT1,KT2
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(4),AL,CP,RA,CK,CE,EPS,PSFC(5)
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSKT2(NZ),BSK(NZ),BSK4(NZ),
     1   BSIT(NX,NY),BSIT2(NX,NY),BSI(NX,NY),BSI4(NX,NY),
     2   BSJT(NX,NY),BSJT2(NX,NY),BSJ(NX,NY),BSJ4(NX,NY)
c
      COMMON/BSAT/ W1(NX,NY,NZ)
      COMMON/BADV/ U1(NX,NY,NZ)
      COMMON/BADV1/ V1(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
c
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
      COMMON/BBA/ XYZ(NX,NY,2),XY1(NX,NY),XY2(NX,NY),Y1(NM),Y2(NM),
     1   Y3(NM),Y4(NM),Y5(NM),Y6(NM),Y10(NM10)
      DIMENSION SQ(1)
      save
C     ******     COMPUTE W1      ***************************************
      A0000=1.
        IF (IMLIFTING .EQ. 0) A0000=0.

      a_irf=1.
        if (irf .eq. 0) a_irf=0.

c$doacross local(j,i)
      DO J=1,JMAX
      DO I=1,IMAX
        XYZ(I,J,1)=0.0
      ENDDO
      ENDDO

      DO 1000 K=2,KLES
        KP=K+1
        RHOK=RHO(K)
        BSKK=BSK(K)
        AMK=AM(K)
        Z3K=Z3(K)
        AM1K=AM1(K)
        RRHO1K=RRHO1(K)
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
         XY1(I,J)=W(I,J,KP)+W(I,J,K)
        XYZ(I,J,2)=RHOK*(XY1(I,J)*XY1(I,J)*RD4Z
     1    -(2.*AK(I,J,K)+BSKK)*(WW1(I,J,KP)-WW1(I,J,K))*AMK*RDZ2
     3    +Z3K*AK(I,J,K)*AK(I,J,K)*RDZ)
          W1(I,J,K)=AM1K*RRHO1K*((XYZ(I,J,1)-XYZ(I,J,2)))
          XYZ(I,J,1)=XYZ(I,J,2)
       ENDDO
       ENDDO

 1000  CONTINUE

       if (ijkadv .eq. 0) then
          y1(2)=0.
          y1(kmax)=0.
          r2dt8=0.008/d2t

          do j=2,jles
          do i=2,iles

            do k=3,kles
               y1(k)=ww1(i,j,k)-ww1(i,j,k-1)
            enddo
            do k=3,kles
               w1(i,j,k)=w1(i,j,k)-r2dt8*(y1(k)-y1(k+1))
            enddo

          enddo
          enddo  !!!sh

       endif  !!sh

cccccc  tao (11-18-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC    *****************************************************************
      IF (IMLIFTING .EQ. 1) THEN
cnew
        DO K=2,KMAX
          Y1K=0.
        DO J=2,JLES
        DO I=2,ILES
          Y1K=Y1K+W1(I,J,K)
        ENDDO
        ENDDO
          Y1K=Y1K*RIJL2
c$doacross local(j,i)       
        DO J=2,JLES
        DO I=2,ILES
           W1(I,J,K)=W1(I,J,K)-Y1K
        ENDDO
        ENDDO
        ENDDO
      ENDIF
C     ******************************************************************
      DO 1500 K=3,KLES
        KP=K+1
        KM=K-1
        KMM=KM-1
C     ***   TURBULENT PROCESSES   **************************************
c$doacross local(j,i,im,y6i)
       DO J=2,JLES
        DO I=2,IMAX
         IM=I-1
         Y6I=AK(I,J,K)+AK(I,J,KM)+AK(IM,J,K)+AK(IM,J,KM)
      XY1(I,J)=-Y6I*((WW1(I,J,K)-WW1(IM,J,K))*R4DX2
     1                   +AM1(K)*(UU1(I,J,K)-UU1(I,J,KM))*RD4Z*RDX)
        ENDDO
       ENDDO

c$doacross local(j,jm,i)
       DO J=2,JMAX
        JM=J-1
        DO I=2,ILES
         Y6(I)=AK(I,J,K)+AK(I,J,KM)+AK(I,JM,K)+AK(I,JM,KM)
         XY2(I,J)=-Y6(I)*((WW1(I,J,K)-WW1(I,JM,K))*R4DY2
     1                   +AM1(K)*(VV1(I,J,K)-VV1(I,J,KM))*RD4Z*RDY)
        ENDDO
       ENDDO

c$doacross local(j,jp,i)
       DO J=2,JLES
        jp=j+1
        DO I=2,ILES
         W1(I,J,K)=W1(I,J,K)+(XY1(I,J)-XY1(I+1,J))+(XY2(I,J)-XY2(I,JP))
        ENDDO
       ENDDO

CCCC     X - COMPONENT     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cnew
c$doacross local(j,i,ipp)
       DO J=2,JLES
        do i=2,iles
          ipp=i+2
          if (i.eq.iles) ipp=3
         xy2(i,j)=a4k(k)*(ww1(ipp,j,k)-3.*(ww1(i+1,j,k)-ww1(i,j,k))
     1                   -ww1(i-1,j,k))
        enddo
        enddo
cnew
c$doacross local(j,i,im)
       DO J=2,JLES
        do i=2,iles
          im=i-1
          if (i.eq.2) im=iles
          w1(i,j,k)=w1(i,j,k)+(xy2(im,j)-xy2(i,j))
        enddo
        enddo

CCCC     Y - COMPONENT     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cnew
c$doacross local(j,jp,jm,jpp,i)
       do j=2,jles
          jp=j+1
          jm=j-1
          jpp=j+2
          if(j.eq.jles) jpp=3
          DO I=2,ILES
         xy1(i,j)=a4k(k)*(ww1(i,jpp,k)-3.*(ww1(i,jp,k)-ww1(i,j,k))
     1                   -ww1(i,jm,k))
        enddo
        enddo
cnew
c$doacross local(j,jm,i)
        do j=2,jles
          jm=j-1
           if (j.eq.2) jm=jles
          DO I=2,ILES
         w1(i,j,k)=w1(i,j,k)+(xy1(i,jm)-xy1(i,j))
        enddo
        enddo
C
       IF (IADVH .EQ. 2) THEN
cnew
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,IMAX
        XYZ(I,J,1)=(U(I,J,K)+U(I,J,KM))*(W(I,J,K)+W(I-1,J,K))*RD4X
       ENDDO
       ENDDO

cnew
c$doacross local(j,jm,i)
       DO J=2,JMAX
       jm=j-1
       DO I=2,ILES
        XYZ(I,J,2)=(V(I,J,K)+V(I,J,KM))*(W(I,J,K)+W(I,JM,K))*RD4Y
       ENDDO
       ENDDO

c$new
c$doacross local(j,jp,i)
       DO J=2,JLES
        jp=j+1
       DO I=2,ILES
         W1(I,J,K)=W1(I,J,K)+(XYZ(I,J,1)-XYZ(I+1,J,1))
     1                     +(XYZ(I,J,2)-XYZ(I,JP,2))
       ENDDO
       ENDDO
       ELSE
CC    ******    4-TH ORDER  HORIZONTAL ADVECTION TERMS   ***************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c$doacross local(j,i,ip,im,ipp)
       do j=2,jles
          do i=2,iles
             ip=i+1
             im=i-1
             ipp=i+2
             if (i.eq.iles) ipp=3
             xy1(i,j)=(u(ip,j,km)+u(ip,j,k))*(d516x*(w(i,j,k)
     1               +w(ip,j,k))-d32x*(w(im,j,k)+w(ipp,j,k)))
          enddo
       enddo

c$doacross local(j,i,im)
       do j=2,jles
        do i=2,imax
          im=i-1
          xy2(i,j)=-d96x*(u(i,j,kmm)+u(i,j,kp))*(w(im,j,k)+w(i,j,k))
        enddo
       enddo

       do j=2,jles
        do i=2,iles
          im=i-1
           if(i.eq.2) im=iles
          ipp=i+2
           if(i.eq.iles) ipp=3
         w1(i,j,k)=w1(i,j,k)+(xy1(im,j)-xy1(i,j)+xy2(im,j)-xy2(ipp,j))
        enddo
       enddo

       DO 160 I=2,ILES
        xy2(i,1)=-d96y*(v(i,1,kmm)+v(i,1,kp))*(w(i,jl2,k)+w(i,1,k))
        xy1(i,1)=(v(i,2,km)+v(i,2,k))*(d516y*(w(i,1,k)+w(i,2,k))
     1                                 -d32y*(w(i,jl2,k)+w(i,3,k)))
        do j=2,jmax
         xy2(i,j)=-d96y*(v(i,j,kmm)+v(i,j,kp))*(w(i,j-1,k)+w(i,j,k))
        enddo
        do j=2,jl2
         xy1(i,j)=(v(i,j+1,km)+v(i,j+1,k))*(d516y*(w(i,j,k)+w(i,j+1,k))
     1                                   -d32y*(w(i,j-1,k)+w(i,j+2,k)))
        enddo
        do j=2,jl2
         w1(i,j,k)=w1(i,j,k)+(xy1(i,j-1)-xy1(i,j)+xy2(i,j-1)-xy2(i,j+2))
        enddo
         w1(i,jles,k)=w1(i,jles,k)+
     1                (xy1(i,jl2)-xy1(i,1)+xy2(i,jl2)-xy2(i,3))
  160  CONTINUE
       ENDIF

Cdan-check
       DO J=2,JLES 
       DO I=2,ILES
         w1(i,j,k)=w1(i,j,k)-a0000*am1(k)*y0(i,j)*
     1                       (wb(kp)*w(i,j,kp)-wb(km)*w(i,j,km))*rd2z
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ENDDO
        ENDDO 
C      
 1500 CONTINUE         

cnew
c$doacross local(j,i)
      do j=2,jles
      do i=2,iles
         w1(i,j,2)=0.
         w1(i,j,kmax)=0.
      enddo
      enddo

C     ******   COMPUTE BUOYANCY AND WATER LOADING TERMS
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
       XY1(I,J)=490.*(DPT(I,J,1)/TB(1)+.61*DQV(I,J,1)+A0000*SQ(1)
     1          -QCL(I,J,1)-QRN(I,J,1)-QCI(I,J,1)-QCS(I,J,1)-QCG(I,J,1))
      ENDDO
      ENDDO

      DO 1801 K=2,KMAX
c$doacross local(j,i,xy2ij)
      DO 1800 J=2,JLES
      DO 180 I=2,ILES
       xy2ij=490.*(DPT(I,J,K)/TB(K)+.61*DQV(I,J,K)+A0000*SQ(K)
     1          -QCL(I,J,K)-QRN(I,J,K)-QCI(I,J,K)-QCS(I,J,K)-QCG(I,J,K))
       W1(I,J,K)=W1(I,J,K)+XY1(I,J)+XY2IJ
       XY1(I,J)=XY2IJ
       W1(I,J,K)=W1(I,J,K)-a_irf*RFA1(K)*WW1(I,J,K) 
  180  continue
 1800  CONTINUE
 1801  CONTINUE
c

      RETURN
      END

      SUBROUTINE ADVECTU (nudge)
C     ****   COMPUTE U1 AND FORCE TERMS
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,ITT=244)
      PARAMETER (NM10=10*NM,itt3=3*itt)
      COMMON /SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
c
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      common/option2/ isngtao,iwbar,iuvbar,isfc,ice,ice2,iradave,idelz
c
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1  NRAN,KT1,KT2
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(4),AL,CP,RA,CK,CE,EPS,PSFC,FCOR,
     1   SEC,AMINUT,RDT
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSKT2(NZ),BSK(NZ),BSK4(NZ),
     1   BSIT(NX,NY),BSIT2(NX,NY),BSI(NX,NY),BSI4(NX,NY),
     2   BSJT(NX,NY),BSJT2(NX,NY),BSJ(NX,NY),BSJ4(NX,NY)
c
      COMMON/BSAT/ W1(NX,NY,NZ)
      COMMON/BADV/ U1(NX,NY,NZ)
      COMMON/BADV1/ V1(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
c
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
c
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
c
      COMMON/BBA/ XYZ(NX,NY,2),XY1(NX,NY),XY2(NX,NY),Y1(NM),Y2(NM),
     1   Y3(NM),Y4(NM),Y5(NM),Y6(NM),Y10(NM10)

      common/bls3/ factor_nuding ! control nudging
      common/bls4a/ ubi(nz),vbi(nz),ub_2h(nz,itt3),vb_2h(nz,itt3)

      real udum1(nz),udum2(nz)
      save
C     ******     COMPUTE U1      ***************************************
C     ******     TOGA-COARE u'w' ***************************************
c      a0000=1.
c        if (nudge .eq. 1) a0000=0.
      a_irf=1.
        if (irf .eq. 0) a_irf=0.
c
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
          XYZ(I,J,1)=0.0
       ENDDO
       ENDDO
C
C-----put the surface momentum fluxes in here
C
       barfct=.5*rho1(2)*100.*100.*rdz  !convert from mks to cgs
       if (itoga .eq. 1) then
          do j=2,jles
             suw(1,j)=suw(iles,j)
          enddo
c$doacross local(j,i)
          do j=2,jles
          do i=2,iles
             xyz(i,j,1)=(suw(i,j)+suw(i-1,j))*barfct
          enddo
          enddo
       endif
CC    ******  
       rdxz4=rd4x*rdz
      DO 2001 K=3,KLES
        KM=K-1
        a11k=AM1(K)*BSK4(K)*R4DZ2
        a22k=AM1(K)*R4DZ2
        a44k=AM(KM)*RRHO(KM)
      DO 2000 J=2,JLES
       DO 210 I=2,ILES
          IM=I-1
          XY1(I,J)=AK(I,J,K)+AK(IM,J,K)+AK(I,J,KM)+AK(IM,J,KM)
         XYZ(I,J,2)=RHO1(K)
     1    *((W(I,J,K)+W(IM,J,K))*(U(I,J,K)+U(I,J,KM))*RD4Z
     2   -a22k*XY1(I,J)*(UU1(I,J,K)-UU1(I,J,KM))
     3   -a11k*(UU1(I,J,K)-UU1(I,J,KM)-UB1(K)+UB1(KM))
     4   -XY1(I,J)*(WW1(I,J,K)-WW1(IM,J,K))*rdxz4)
        U1(I,J,KM)=a44k*(XYZ(I,J,1)-XYZ(I,J,2))
         XYZ(I,J,1)=XYZ(I,J,2)
  210   CONTINUE
 2000  CONTINUE
 2001  CONTINUE

       A0R=AM(KLES)*RRHO(KLES)
cc$doacross local(j,i,a0r)
       DO J=2,JLES
       DO I=2,ILES
          U1(I,J,KLES)=A0R*XYZ(I,J,2)
       ENDDO
       ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-18-97)
c !!shie 10-14-98
       if (ijkadv .eq. 0) then

          y1(2)=0.
          y1(kmax)=0.

          r2dt8=0.008/d2t

          do j=2,jles
          do i=2,iles

            do k=3,kles
               y1(k)=uu1(i,j,k)-uu1(i,j,k-1)-ub1(k)+ub1(k-1)
            enddo
            do k=2,kles
               u1(i,j,k)=u1(i,j,k)-r2dt8*(y1(k)-y1(k+1))
            enddo

          enddo
          enddo
        endif

cccccc  tao (11-18-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (nudge .ne. 1) THEN
        DO K=2,Kmax
          Y1(K)=0.
        ENDDO
c$doacross local(k,j,i)
        DO K=2,KLES
        DO J=2,JLES
        DO I=2,ILES
         Y1(K)=Y1(K)+U1(I,J,K)
        ENDDO
        ENDDO
        ENDDO
        DO K=2,KLES
          Y1(K)=Y1(K)*RIJL2
        ENDDO
cnew
c$doacross local(k,j,i)
        DO K=2,KLES
        DO J=2,JLES
        DO I=2,ILES
          U1(I,J,K)=U1(I,J,K)-Y1(K)
        ENDDO
        ENDDO
        ENDDO
       ENDIF
CC    ******              HORIZONTAL ADVECTION TERMS   *****************
       a11k=RDX*RD4Y

      DO 2500 K=2,KLES
       KM=K-1
        z33k=Z3(K)*RDX
        a22k=AM(K)*RD4Z
           ub2=0.
           vb2=0.
C     ****    TURBULENT PROCESSES    ***********************************
c$doacross local(j,i)
       DO J=2,JLES
        DO I=1,ILES
          XY1(I,J)=-2.*AK(I,J,K)*(UU1(I+1,J,K)-UU1(I,J,K))*RDX2
     1            +Z33K*AK(I,J,K)*AK(I,J,K)
        ENDDO
       ENDDO

c$doacross local(j,jp,i)
       DO J=1,JLES
        jp=j+1
        DO I=2,ILES
          Y6(I)=AK(I,J,K)+AK(I-1,J,K)+AK(I,JP,K)+AK(I-1,JP,K)
         XY2(I,J)=-Y6(I)*((UU1(I,JP,K)-UU1(I,J,K))*R4DY2
     1            +(VV1(I,JP,K)-VV1(I-1,JP,K))*a11k)
        ENDDO
       ENDDO

c$doacross local(j,jm,i)
        DO J=2,JLES
           jm=j-1
        DO I=2,ILES
         U1(I,J,K)=U1(I,J,K)+(XY1(I-1,J)-XY1(I,J))
     1                      +(XY2(I,JM)-XY2(I,J))
        ENDDO
        ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$doacross local(j,i,ipp)
       DO J=2,JLES
        do i=2,iles
         ipp=i+2
          if (i .eq. iles) ipp=3
         xy1(i,j)=a4k(k)*(uu1(ipp,j,k)-3.*(uu1(i+1,j,k)-uu1(i,j,k))
     1                   -uu1(i-1,j,k))
        enddo
        enddo
c$doacross local(j,i,im)
       DO 840 J=2,JLES
        do i=2,iles
          im=i-1
           if (i .eq. 2) im=iles
         u1(i,j,k)=u1(i,j,k)+(xy1(im,j)-xy1(i,j))
       enddo
  840  CONTINUE

c$doacross local(j,jm,jp,jpp,i)
        do j=2,jles
          jpp=j+2
          if (j .eq. jles) jpp=3
          jp=j+1
          jm=j-1
         DO I=2,ILES
        xy1(i,j)=a4k(k)*(uu1(i,jpp,k)-3.*(uu1(i,jp,k)-uu1(i,j,k))
     1                   -uu1(i,jm,k))
        enddo
        enddo

c$doacross local(j,jm,i)
       do j=2,jles
          jm=j-1
          if(j.eq.2) jm=jles
       DO I=2,ILES
          u1(i,j,k)=u1(i,j,k)+(xy1(i,jm)-xy1(i,j))
         enddo
         enddo
C
       IF (IADVH .EQ. 2) THEN
       DO J=1,JLES
       DO I=2,ILES
       XYZ(I,J,2)=(U(I,J,K)+U(I,J+1,K))
     1             *(V(I,J+1,K)+V(I-1,J+1,K)-VB2)*RD4Y
       ENDDO
       ENDDO
       DO J=2,JLES
       DO I=1,ILES
         XYZ(I,J,1)=(U(I+1,J,K)+U(I,J,K))*(U(I+1,J,K)+U(I,J,K)-UB2)*RD4X
       ENDDO
       ENDDO
       DO J=2,JLES
       DO I=2,ILES
        U1(I,J,K)=U1(I,J,K)+(XYZ(I-1,J,1)-XYZ(I,J,1))
     1                     +(XYZ(I,J-1,2)-XYZ(I,J,2))
       ENDDO
       ENDDO
       ELSE
C     ******    4-TH ORDER ADV SCHEME      *****************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$doacross local(j,i,ipp,ip,im)
       DO J=2,JLES
        do i=2,iles
          ipp=i+2
           if (i.eq.iles) ipp=3
          ip=i+1
          im=i-1
         xy1(i,j)=d516x*(u(i,j,k)+u(ip,j,k)-ub2)*(u(i,j,k)+u(ip,j,k))
     1           -d32x*(u(im,j,k)+u(ipp,j,k)-ub2)
     2                 *(u(im,j,k)+u(ipp,j,k))
        enddo
        enddo

c$doacross local(j,i,im)
       DO J=2,JLES
        do i=2,imax
         im=i-1
         xy2(i,j)=-d96x*(u(im,j,k)+u(i,j,k)-ub2)*(u(im,j,k)+u(i,j,k))
        enddo
        enddo

c$doacross local(j,i,im,ipp)
       DO 240 J=2,JLES
        do i=2,iles
          im=i-1
          if (i.eq.2) im=iles
          ipp=i+2
          if (i.eq.iles) ipp=3
         u1(i,j,k)=u1(i,j,k)+(xy1(im,j)-xy1(i,j)+xy2(im,j)-xy2(ipp,j))
        enddo
  240  CONTINUE
C

c$doacross local(j,jm,jp,jpp,i,imm)
       do j=2,jles
          jm=j-1
          jp=j+1
          jpp=j+2
          if (j.eq.jles) jpp=3
          DO I=2,ILES
             imm=i-2
             if (i.eq.2) imm=il2
             xy1(i,j)=(v(i-1,jp,k)+v(i,jp,k)-vb2)*
     1                  (d516y*(u(i,j,k)+u(i,jp,k))
     1                 -d32y*(u(i,jm,k)+u(i,jpp,k)))
        enddo
        enddo

c$doacross local(j,jm,i,imm)
      do j=2,jmax
        jm=j-1
        DO I=2,ILES
           imm=i-2
           if (i.eq.2) imm=il2
         xy2(i,j)=-d96y*(u(i,jm,k)+u(i,j,k))*(v(i+1,j,k)+v(imm,j,k)-vb2)
        enddo
      enddo

c$doacross local(j,jpp,jm,i)
      do j=2,jles
         jpp=j+2
         if (j .eq. jles) jpp=3
         jm=j-1
         if(j.eq.2) jm=jles
         DO I=2,ILES
           u1(i,j,k)=u1(i,j,k)+(xy1(i,jm)-xy1(i,j)+xy2(i,jm)-xy2(i,jpp))
       enddo
       enddo
       ENDIF

c$doacross local(j,i)
       DO 881 J=2,JLES
       DO 880 I=2,ILES
          U1(I,J,K)=U1(I,J,K)-a_irf*rfa(k)*(uu1(i,j,k)-ub1(k))
  880  CONTINUE
  881  CONTINUE
 2500  CONTINUE

c
c
cshie 9/12/98
      if(nudge.eq.1) then
         do k=2,kles
         do j=2,jles
         do i=2,iles
           u1(i,j,k)=u1(i,j,k)-factor_nuding*(ub(k)-ubi(k)) !!sh as in 2D
         enddo
         enddo
         enddo
      endif
cshie 9/12/98
c
      if (isec.eq.isec/86400*86400) then

         print *,'!!! for u(i,k),ubi(k),isec=',isec
         do k=2,kles
          udum1(k)=0.
          udum2(k)=0.
         do j=2,jles
         do i=2,iles
          udum1(k)=udum1(k)+u1(i,j,k)*30.+u(i,j,k)
          udum2(k)=udum2(k)+u(i,j,k)
         enddo
         enddo
         enddo

         do k=2,kles
          udum1(k)=udum1(k)/(real(iles-1)*real(jles-1))
          udum2(k)=udum2(k)/(real(iles-1)*real(jles-1))
         enddo

      WRITE(6,2034)
      DO K1=2,KMAX
        K2=KMAX+2-K1
      WRITE(6,2031) K2,udum2(K2),udum1(K2),UBi(K2)
        enddo

        endif

 2031 FORMAT(1X,I4,3F10.2)
 2034 FORMAT(//,'  LEVEL  UBm1    UBp1   UBi')

      if (nudge .ne. 1) then
c$doacross local(k,ubtk,j,i)
        do 4002 k=2,kles
         ubtk=ubt(k)
        do 4001 j=2,jles
        do 4000 i=2,iles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           u1(i,j,k)=u1(i,j,k)+ubtk
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 4000   continue
 4001   continue
 4002   continue
      endif
      RETURN
      END

      SUBROUTINE ADVECTV (nudge)
C     ****   COMPUTE V1 AND FORCE TERMS
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,ITT=244)
      PARAMETER (NM10=10*NM,itt3=3*itt)
      COMMON /SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
c
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      common/option2/ isngtao,iwbar,iuvbar,isfc,ice,ice2,iradave,idelz
c
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1  NRAN,KT1,KT2
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(4),AL,CP,RA,CK,CE,EPS,PSFC,FCOR,
     1   SEC,AMINUT,RDT
c
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSKT2(NZ),BSK(NZ),BSK4(NZ),
     1   BSIT(NX,NY),BSIT2(NX,NY),BSI(NX,NY),BSI4(NX,NY),
     2   BSJT(NX,NY),BSJT2(NX,NY),BSJ(NX,NY),BSJ4(NX,NY)
c
      COMMON/BSAT/ W1(NX,NY,NZ)
      COMMON/BADV/ U1(NX,NY,NZ)
      COMMON/BADV1/ V1(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
c
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/BLS/ Y0(NX,NY),TS0(NX,NY),QSS0(NX,NY)
c
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
c
      COMMON/BBA/ XYZ(NX,NY,2),XY1(NX,NY),XY2(NX,NY),Y1(NM),Y2(NM),
     1   Y3(NM),Y4(NM),Y5(NM),Y6(NM),Y10(NM10)

      common/bls3/ factor_nuding ! control nudging
      common/bls4a/ ubi(nz),vbi(nz),ub_2h(nz,itt3),vb_2h(nz,itt3)
      save
C     ******      COMPUTE V1      **************************************
C     ******      TOGA-COARE v'w' **************************************      
c      A0000=1.
c        IF (nudge .EQ. 1) A0000=0.
      a_irf=1.
        if (irf .eq. 0) a_irf=0.
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
          XYZ(I,J,1)=0.0
       ENDDO
       ENDDO
C
C-----put the surface momentum fluxes in here
C
       BARFCT=.5*RHO1(2)*100.*100.*RDZ  !convert from mks to cgs
        IF (ITOGA .EQ. 1) THEN
          DO J=2,JLES
             JM=J-1
              IF(J.EQ.2) JM=JLES
          DO I=2,ILES
             XYZ(I,J,1)=(SVW(I,J)+SVW(I,JM))*BARFCT
          ENDDO
          ENDDO
        endif
CC    ******        
       rdyz4=RDY*RD4Z
      DO 3001 K=3,KLES
        KM=K-1
        a11k=AM1(K)*R4DZ2
        a22k=AM1(K)*BSK4(K)*R4DZ2
        a44k=AM(KM)*RRHO(KM)
      DO 3000 J=2,JLES
       DO 310 I=2,ILES
        XY1(I,J)=AK(I,J,K)+AK(I,J,KM)+AK(I,J-1,K)+AK(I,J-1,KM)
        XYZ(I,J,2)=RHO1(K)*((W(I,J,K)+W(I,J-1,K))
     1                      *(V(I,J,KM)+V(I,J,K))*RD4Z
     2                     -a11k*XY1(I,J)*(VV1(I,J,K)-VV1(I,J,KM))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     3                     -a22k*(VV1(I,J,K)-VV1(I,J,KM)-VB1(K)+VB1(KM))
     4                     -XY1(I,J)*(WW1(I,J,K)-WW1(I,J-1,K))*rdyz4)
        V1(I,J,KM)=a44k*(XYZ(I,J,1)-XYZ(I,J,2))
         XYZ(I,J,1)=XYZ(I,J,2)
  310  CONTINUE
 3000 CONTINUE
 3001 CONTINUE
        A0R=AM(KLES)*RRHO(KLES)
c$doacross local(j,i)
      DO 321 J=2,JLES
      DO 320 I=2,ILES
       V1(I,J,KLES)=A0R*XYZ(I,J,2)
  320  CONTINUE
  321  CONTINUE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-18-97)
      if (ijkadv .eq. 0) then
          y2(2)=0.
          y2(kmax)=0.
          r2dt8=0.008/d2t
          do j=2,jles
          do i=2,iles
            do k=3,kles
             y2(k)=vv1(i,j,k)-vv1(i,j,k-1)-vb1(k)+vb1(k-1)
            enddo
            do k=2,kles
             v1(i,j,k)=v1(i,j,k)-r2dt8*(y2(k)-y2(k+1))
            enddo
          enddo
          enddo
      endif
cccccc  tao (11-18-97)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if (nudge .ne. 1) then
        DO K=2,Kmax
          Y1(K)=0.
        ENDDO
c$doacross local(k,j,i)
        DO K=2,KLES
        DO J=2,JLES
        DO I=2,ILES
         Y1(K)=Y1(K)+V1(I,J,K)
        ENDDO
        ENDDO
        ENDDO
        DO K=2,KLES
          Y1(K)=Y1(K)*RIJL2
c$doacross local(j,i)
        DO J=2,JLES
        DO I=2,ILES
          V1(I,J,K)=V1(I,J,K)-Y1(K)
        ENDDO
        ENDDO
        ENDDO
       endif
CC    ******     HORIZONTAL ADVECTION TERMS   ***************************
       rrdxy4=RDX*RD4Y
      DO 3500 K=2,KLES
        KM=K-1
        z33k=Z3(K)*RDY
        a11k=AM(K)*RD4Z
          ub2=0.
          VB2=0.
C     ****   TURBULENT PROCESSES        ********************************
c$doacross local(j,jm,i,ip,y6i)
       DO J=2,JLES
        jm=j-1
        DO I=1,ILES
         ip=i+1
         Y6I=AK(I,J,K)+AK(IP,J,K)+AK(I,JM,K)+AK(IP,JM,K)
        XY1(I,J)=-Y6I*((VV1(IP,J,K)-VV1(I,J,K))*R4DX2
     1                  +(UU1(IP,J,K)-UU1(IP,JM,K))*rrdxy4)
        ENDDO
       ENDDO

c$doacross local(j,i)
       DO J=1,JLES
        DO I=2,ILES
         XY2(I,J)=-2.*AK(I,J,K)*(VV1(I,J+1,K)-VV1(I,J,K))*RDY2
     1            +z33k*AK(I,J,K)*AK(I,J,K)
        ENDDO
       ENDDO

c$doacross local(j,jm,i)
        DO J=2,JLES
          JM=J-1
        DO I=2,ILES
         V1(I,J,K)=V1(I,J,K)+(XY1(I-1,J)-XY1(I,J))+(XY2(I,JM)-XY2(I,J))
        ENDDO
        ENDDO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$doacross local(j,i,ipp)
       DO J=2,JLES
        do i=2,iles
          ipp=i+2
           if (i.eq.iles) ipp=3
        xy1(i,j)=a4k(k)*(vv1(ipp,j,k)-3.*(vv1(i+1,j,k)-vv1(i,j,k))
     1                   -vv1(i-1,j,k))
        enddo
        enddo
c$doacross local(j,i,im)
       DO 520 J=2,JLES
        do i=2,iles
          im=i-1
           if (i.eq.2) im=iles
        v1(i,j,k)=v1(i,j,k)+(xy1(im,j)-xy1(i,j))
        enddo
  520  CONTINUE

c$doacross local(j,jpp,jp,jm,i)
       do j=2,jles
          jpp=j+2
          if (j.eq.jles) jpp=3
          jp=j+1
          jm=j-1
          DO I=2,ILES
             xy1(i,j)=a4k(k)*(vv1(i,jpp,k)-3.*(vv1(i,j+1,k)-vv1(i,j,k))
     1               -vv1(i,j-1,k))
        enddo
        enddo

c$doacross local(j,jm,i)
       do j=2,jles
          jm=j-1
          if(j.eq.2) jm=jles
          DO I=2,ILES
             v1(i,j,k)=v1(i,j,k)+(xy1(i,jm)-xy1(i,j))
          enddo
       enddo

       IF (IADVH .EQ. 2) THEN
       DO J=2,JLES
         JM=J-1
       DO I=1,ILES
        XYZ(I,J,1)=(V(I,J,K)+V(I+1,J,K))*
     1             (U(I+1,J,K)+U(I+1,JM,K)-UB2)*RD4X
       ENDDO
       ENDDO

c$doacross local(j,jp,i)
       DO J=1,JLES
        JP=J+1
       DO I=2,ILES
        XYZ(I,J,2)=(V(I,J,K)+V(I,JP,K))*(V(I,J,K)+V(I,JP,K)-VB2)*RD4Y
       ENDDO
       ENDDO

c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
         V1(I,J,K)=V1(I,J,K)+(XYZ(I-1,J,1)-XYZ(I,J,1))
     1                     +(XYZ(I,J-1,2)-XYZ(I,J,2))
       ENDDO
       ENDDO
       ELSE
CCCC  ******    4-TH ORDER HORIZONTAL ADV SCHEME    ********************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$doacross local(j,jmm,jm,i,ipp,ip)
       DO J=2,JLES
         JMM=J-2
          IF (J.EQ.2) JMM=JL2
          jm=j-1
          do i=2,iles
             ipp=i+2
             if (i.eq.iles) ipp=3
             ip=i+1
             xy1(i,j)=(u(ip,j,k)+u(ip,jm,k)-ub2)
     1               *(d516x*(v(i,j,k)+v(ip,j,k))
     1               -d32x*(v(i-1,j,k)+v(ipp,j,k)))
        enddo
        enddo

c$doacross local(j,jmm,jp,i)
      DO J=2,JLES
         JMM=J-2
         IF (J.EQ.2) JMM=JL2
         jp=j+1
         do i=2,imax
            xy2(i,j)=-d96x*(u(i,jp,k)+u(i,jmm,k)-ub2)*
     1                  (v(i,j,k)+v(i-1,j,k))
         enddo
      enddo

c$doacross local(j,i,im,ipp)
       DO J=2,JLES
          do i=2,iles
             im=i-1
             if (i.eq.2) im=iles
             ipp=i+2
             if(i.eq.iles) ipp=3
         v1(i,j,k)=v1(i,j,k)+(xy1(im,j)-xy1(i,j)+xy2(im,j)-xy2(ipp,j))
       enddo
       enddo

c$doacross local(j,jm,jp,jpp,i)
       do j=1,jles
          jm=j-1
          if(j.eq.1) jm=jl2
          jp=j+1
          jpp=j+2
          if(j.eq.jles) jpp=3
          DO I=2,ILES
            xy1(i,j)=d516y*(v(i,j,k)+v(i,jp,k)-vb2)*(v(i,j,k)+v(i,jp,k))
     1          -d32y*(v(i,jm,k)+v(i,jpp,k)-vb2)*(v(i,jm,k)+v(i,jpp,k))
          enddo
       enddo

c$doacross local(j,jm,i)
      do j=1,jmax
         jm=j-1
         if(j.eq.1) jm=jl2
         DO I=2,ILES
         xy2(i,j)=-d96y*(v(i,jm,k)+v(i,j,k)-vb2)*(v(i,jm,k)+v(i,j,k))
         enddo
      enddo

c$doacross local(j,jpp,jm,i)
      do j=2,jles
         jpp=j+2
         if (j.eq.jles) jpp=3
         jm=j-1
         if (j.eq.2) jm=jles
         DO I=2,ILES
           v1(i,j,k)=v1(i,j,k)+(xy1(i,jm)-xy1(i,j)+xy2(i,jm)-xy2(i,jpp))
         enddo
      enddo

       ENDIF

c$doacross local(j,i)
       DO 561 J=2,JLES
       DO 560 I=2,ILES
           V1(I,J,K)=V1(I,J,K)-a_irf*rfa(k)*(vv1(i,j,k)-vb1(k))
  560  CONTINUE
  561  CONTINUE

 3500  CONTINUE

cshie 9/12/98
      if (nudge .eq. 1) then
         do k=2,kles
         do j=2,jles
         do i=2,iles
           v1(i,j,k)=v1(i,j,k)-factor_nuding*(vb(k)-vbi(k)) !!sh as in 2D
         enddo
         enddo
         enddo
      endif
cshie 9/12/98

      if (nudge .ne. 1) then
c$doacross local(k,vbtk,j,i)
        do 4002 k=2,kles
        vbtk=vbt(k)
        do 4001 j=2,jles
        do 4000 i=2,iles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           v1(i,j,k)=v1(i,j,k)+vbtk
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 4000   continue
 4001   continue
 4002   continue
      endif
c$doacross local(k,j,i)
      DO K=1,KMAX
      DO J=1,JMAX
      DO I=1,IMAX
       AK(I,J,K)=2.*AK(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END


      SUBROUTINE CMPF
C   ****   COMPUTE FORCE TERMS
C   ****   UPDATE W, U AND V WITHOUT PRESSURE FORCE
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      parameter (nm14=14*nm)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
c
      COMMON/BSAT/ W1(NX,NY,NZ)
      COMMON/BADV/ F(NX,NY,NZ)
      COMMON/BADV1/ V1(NX,NY,NZ)
      COMMON/B3Uv/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B3Wp/ W(NX,NY,NZ)
      COMMON/B4Uv/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4Wp/ WW1(NX,NY,NZ)
c
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      common/b555/ BA(NZ),BB(NZ)
      COMMON/BBA/ XY1(NX,NY),XY2(NX,NY),XY3(NX,NY),XY4(NX,NY),Y1(NM),
     1   Y2(NM),y3(nm14)
      common/TMP/ y1d(nx,ny,nz),y2d(nx,ny,nz),y3d(nx,ny,nz),
     1            y4d(nx,ny,nz),y03d(nx,ny,1)
c      save rdt22
      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rdt22=1./d2t
      call bndop (f,v1)

      IF (LIPPS .eq. 0) THEN

      DO K=2,KLES
         A1=AM(K)*RDZ*RRHO(K)
         A2=DZ2/(CP*TB(K))
        DO J=2,JLES
          DO I=2,ILES
c           Y2(I)=F(I,J,K)
            Y4D(I,J,K)=F(I,J,K)
            F(I,J,K)=A2*(( F(i+1,j,k) - F(i,j,k))*RDX
     1         +(V1(I,J+1,K)-V1(I,J,K))*RDY
     2         +A1*(W1(I,J,K+1)*RHO1(K+1)-W1(I,J,K)*RHO1(K))
     3         +((UU1(i+1,J,K)-UU1(I,J,K))*RDX
     4         +(VV1(I,J+1,K)-VV1(I,J,K))*RDY
     5         +(WW1(I,J,K+1)*RHO1(K+1)-WW1(I,J,K)*RHO1(K))*A1)*rdt22)
c            Y1(I)=U(I,J,K)+EPS*(-2.*U(I,J,K)+UU1(I,J,K))
             Y3D(I,J,K)=U(I,J,K)+EPS*(-2.*U(I,J,K)+UU1(I,J,K))
c            U(I,J,K)=UU1(I,J,K)+Y2(I)*D2T
c            UU1(I,J,K)=Y1(I)
          ENDDO
          ENDDO
          ENDDO

      ELSE

c$doacross local(k,kp,a1,a2,j,jp,i,ip,yy2)
      DO K=2,KLES
         KP=K+1
         A1=AM(K)*RDZ*RRHO(K)
         a2=dz2/cp
         DO J=2,JLES
            JP=J+1
            DO I=2,ILES
               IP=I+1
               Y4D(I,J,K)=F(I,J,K)
               F(I,J,K)=A2*((F(ip,j,k)-F(i,j,k))*RDX
     1                 +(V1(I,JP,K)-V1(I,J,K))*RDY
     2                 +A1*(W1(I,J,KP)*RHO1(KP)-W1(I,J,K)*RHO1(K))
     3                 +((UU1(IP,J,K)-UU1(I,J,K))*RDX
     4                 +(VV1(I,JP,K)-VV1(I,J,K))*RDY
     5             +(WW1(I,J,KP)*RHO1(KP)-WW1(I,J,K)*RHO1(K))*A1)*rdt22)
               Y3D(I,J,K)=U(I,J,K)+EPS*(-2.*U(I,J,K)+UU1(I,J,K))
            ENDDO
         ENDDO
      ENDDO

      ENDIF

c$doacross local(k,j,i)
      DO K=3,KLES
      DO J=2,JLES
        DO I=2,ILES
             y1d(i,j,k)=W(I,J,K)+EPS*(-2.*W(I,J,K)+WW1(I,J,K))
           W(I,J,K)=WW1(I,J,K)+W1(I,J,K)*D2T
         ENDDO
        ENDDO
      ENDDO

c$doacross local(k,j,i)
      DO K=3,KLES
      DO J=2,JLES
        DO I=2,ILES
           WW1(I,J,K)=Y1D(I,J,K)
         ENDDO
        ENDDO
      ENDDO

c$doacross local(k,j,i)
      DO K=2,KLES
      DO J=2,JLES
        DO I=2,ILES
             Y2D(I,J,K)=V(I,J,K)+EPS*(-2.*V(I,J,K)+VV1(I,J,K))
           V(I,J,K)=VV1(I,J,K)+V1(I,J,K)*D2T
           U(I,J,K)=UU1(I,J,K)+Y4D(I,J,K)*D2T
         ENDDO
        ENDDO
      ENDDO

c$doacross local(k,j,i)
      DO K=2,KLES
      DO J=2,JLES
        DO I=2,ILES
           VV1(I,J,K)=Y2D(I,J,K)
           UU1(I,J,K)=Y3D(I,J,K)
         ENDDO
        ENDDO
      ENDDO

C   ****   SET BOUNDARY CONDITIONS
      IF (LIPPS .EQ. 1) THEN
        A1=DZ*BA(2)/(CP*AM1(2))
        A2=DZ*BB(KLES)/(CP*AM1(KLES))
      ELSE
        A1=DZ*BA(2)/(CP*TB(2)*AM1(2))
        A2=DZ*BB(KLES)/(CP*TB(KLES)*AM1(KLES))
      ENDIF

c$doacross local(j,i)
      DO 400 J=2,JLES
      DO 40 I=2,ILES
        F(I,J,2)=F(I,J,2)+W1(I,J,2)*A1
        F(I,J,KLES)=F(I,J,KLES)-W1(I,J,KMAX)*A2
   40 CONTINUE
  400 CONTINUE

      RETURN
      END


      SUBROUTINE CMPWUV
C   ****   COMPUTE NEW W, U AND V
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX)
      PARAMETER (ND1=6*NM)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
c
      COMMON/BADV/ PI(NX,NY,NZ)
      COMMON/B3Uv/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B3Wp/ W(NX,NY,NZ)
      COMMON/B4Uv/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4Wp/ WW1(NX,NY,NZ)
c
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1   QA1(NZ),COEF(NZ),C1(NZ),C2(NZ),C3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2   VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/BA/ Y1(NM),Y2(NM),Y3(NM),Y4(NM),XXW(ND1)
      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF (LIPPS .EQ. 1) THEN
         A2=CP*D2T*RDZ
c$doacross local(k,km,yy1,j,i)
         DO K=3,KLES
            KM=K-1
            YY1=AM1(K)*A2
            DO J=2,JLES
               DO I=2,ILES
                  W(I,J,K)=W(I,J,K)-YY1*(PI(I,J,K)-PI(I,J,KM))
                  WW1(I,J,K)=WW1(I,J,K)+EPS*W(I,J,K)
               ENDDO
            ENDDO
          ENDDO
 
      ELSE

         A2=.5*CP*D2T*RDZ
         DO K=3,KLES
            KM=K-1
            Y1(K)=AM1(K)*(TB(KM)+TB(K))*A2
            DO J=2,JLES
               DO I=2,ILES
                  W(I,J,K)=W(I,J,K)-Y1(K)*(PI(I,J,K)-PI(I,J,KM))
                  WW1(I,J,K)=WW1(I,J,K)+EPS*W(I,J,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      CALL BNDOP (PI,PI)
      A1=CP*D2T*RDX
      A2=CP*D2T*RDY
      DO 15 K=2,KLES
       IF (LIPPS .EQ. 1) THEN
        Y1(K)=A1
        Y2(K)=A2
       ELSE
        Y1(K)=TB(K)*A1
        Y2(K)=TB(K)*A2
       ENDIF
   15 CONTINUE
c$doacross local(k,j,i)
      DO 301 K=2,KLES
      DO 300 J=2,JLES
       DO 30 I=2,ILES
         U(I,J,K)=U(I,J,K)-Y1(K)*(PI(I,J,K)-PI(I-1,J,K))
         UU1(I,J,K)=UU1(I,J,K)+EPS*U(I,J,K)
         V(I,J,K)=V(I,J,K)-Y2(K)*(PI(I,J,K)-PI(I,J-1,K))
         VV1(I,J,K)=VV1(I,J,K)+EPS*V(I,J,K)
   30  CONTINUE
  300 CONTINUE
  301 CONTINUE
c$doacross local(j,i)
      DO 401 J=2,JLES
      DO 400 I=2,ILES
       W(I,J,2)=0.
       WW1(I,J,2)=0.
  400 CONTINUE
  401 CONTINUE
      CALL bndop (U,UU1)
      CALL bndop (V,VV1)
      CALL bndop (W,WW1)
c$doacross local(j,i)
       DO J=1,JMAX
       DO I=1,IMAX
        W(I,J,2)=0.
        WW1(I,J,2)=0.
        W(I,J,KMAX)=0.0
        WW1(I,J,KMAX)=0.0
       ENDDO
       ENDDO

      DO K=1,KMAX
       Y1(K)=0.
       Y2(K)=0.
       Y3(K)=0.
       Y4(K)=0.
      ENDDO

c$doacross local(k,j,i)
      DO K=1,KMAX
      DO 511 J=2,JLES
      DO 510 I=2,ILES
        Y1(K)=Y1(K)+U(I,J,K)
        Y2(K)=Y2(K)+UU1(I,J,K)
        Y3(K)=Y3(K)+V(I,J,K)
        Y4(K)=Y4(K)+VV1(I,J,K)
  510 CONTINUE
  511 CONTINUE
       ENDDO

      DO K=1,KMAX
        UB(K)=Y1(K)*RIJL2
        UB1(K)=Y2(K)*RIJL2
        VB(K)=Y3(K)*RIJL2
        VB1(K)=Y4(K)*RIJL2
       ENDDO
      RETURN
      END

      SUBROUTINE SLVPI (IFLG)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****  SOLVE 3-D PRESSURE EQUATION FOR CYCLIC BOUNDARY CONDITIONS C
C     ****   STARTS ON APRIL 1996                                      C
C     ****   CALLS ROUTINES FFTX & FFTY                                C
C     ****   CALLS SLVFY & SLVBY IF NX DOES NOT EQUAL NY               C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NXI=NX-2,NYJ=NY-2,NZK=NZ-2)
      PARAMETER (ND1=NX*NY*NZ-NXI*NYJ*NZK-2*NX*NZ-NX)
      PARAMETER (ND2=NX*NY*NZ-NXI*NYJ*NZK-NX*NZ)
      COMMON/RTIME/ ISEC1
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
c
      COMMON/BSAT/ AR(NZK,NYJ,NXI),EE(NX,NZ),FF(NX,NZ),Y1(NX),XXD(ND1)
      COMMON/BSAT1/ AI(NZK,NYJ,NXI),AUX(NX,NZ),YYD(ND2)
      COMMON/BADV/ F(NX,NY,NZ)
      COMMON/B555/ CC(NZ),AA(NZ)
c
      COMMON/BSLVPI/ BETA(NX),BETA1(NY),BT(NXI,NYJ)
      save
      IF(IFLG.EQ.1) GO TO 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C2PI=8.*ATAN(1.)
      RIL2=1./FLOAT(IL2)
      IHF=IL2/2+1
      IHF1=IHF+1
      IHF2=IHF1+1
      PIX=C2PI*RIL2
      DO I=2,ILES
        BETA(I)=2.*DZ2*(COS(PIX*FLOAT(I-2))-1.)*RDX2
      ENDDO
       BETA(1)=0.
       BETA(IMAX)=0.
C     ******      ******
      CALL FFTX (0,IL2)
C     ******      ******
      RJL2=1./FLOAT(JL2)
      JHF=JL2/2+1
      JHF1=JHF+1
      JHF2=JHF1+1
      PIY=C2PI*RJL2
      DO J=2,JLES
         BETA1(J)=2.*DZ2*(COS(PIY*FLOAT(J-2))-1.)*RDY2
      ENDDO
       BETA1(1)=0.
       BETA1(JMAX)=0.
C     ******      ******
      CALL FFTY (0,JL2)
C     ******      ******

c$doacross local(j,jm,i,im)
      DO J=2,JLES
      JM=J-1
      DO I=2,ILES
        IM=I-1
        BT(IM,JM)=BETA(I)+BETA1(J)
      ENDDO
      ENDDO
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ***    FORWARD TRANSFORMATION IN X AND Y    **********************
C     ***    FORWARD TRANSFORMATION IN X

    1 continue
c$doacross local(I,IM,J,JM,K)
      DO I=2,ILES
      IM=I-1
      DO J=2,JLES
      JM=J-1
      DO K=2,KLES
       AR(K-1,JM,IM)=F(I,J,K)*RIL2
      ENDDO
      ENDDO
      ENDDO

c$doacross local(i,j,k)
      DO I=1,IL2
      DO J=1,JL2
      DO K=1,KL2
        AI(K,J,I)=0.0
      ENDDO
      ENDDO
      ENDDO

C     ******      ******
      CALL FFTX (1,IL2)
C     ******      ******

c$doacross local (k,km,j,jm,i)
      DO K=2,KLES
      KM=K-1
      DO J=2,JLES
        JM=J-1
        DO I=2,IHF1
         F(I,J,K)=AR(KM,JM,I-1)
        ENDDO
      ENDDO
      ENDDO

c$doacross local (k,km,j,jm,i)
      DO K=2,KLES
      KM=K-1
      DO J=2,JLES
        JM=J-1
        DO I=IHF2,ILES
         F(I,J,K)=AI(KM,JM,I-1)
        ENDDO
      ENDDO
      ENDDO

C     ***    FORWARD TRANSFORMATION IN Y    ***************************
      IF (IL2 .EQ. JL2) THEN

c$doacross local(J,JM,I,IM,K)
        DO J=2,JLES
        JM=J-1
        DO I=2,ILES
        IM=I-1
        DO K=2,KLES
          AR(K-1,IM,JM)=F(I,J,K)*RJL2
        ENDDO
        ENDDO
        ENDDO
c$doacross local(j,i,k)
        DO J=1,JL2
        DO I=1,IL2
        DO K=1,KL2
          AI(K,I,J)=0.0
        ENDDO
        ENDDO
        ENDDO
C     ******      ******
       CALL FFTY (1,JL2)
C     ******      ******

c$doacross local(k,km,j,jm,i)
        DO K=2,KLES
          KM=K-1
          DO J=2,JHF1
            JM=J-1
            DO I=2,ILES
              F(I,J,K)=AR(KM,I-1,JM)
            ENDDO
          ENDDO
        ENDDO

c$doacross local(k,km,j,jm,i)
        DO K=2,KLES
          KM=K-1
            DO J=JHF2,JLES
              JM=J-1
              DO I=2,ILES
                F(I,J,K)=AI(KM,I-1,JM)
            ENDDO
          ENDDO
        ENDDO

      ELSE
C     ***    FORWARD TRANSFORMATION IN Y        ************************
        CALL SLVFY (JHF1,JHF2,RJL2)
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ******************************************************************
      DO 3000 J=2,JLES
        DO I=2,ILES
         AUX(I,2)=1./(AA(2)-BT(I-1,J-1))
         EE(I,2)=AA(2)*AUX(I,2)
         FF(I,2)=-F(I,J,2)*AUX(I,2)
        ENDDO
       DO 320 K=3,KL2
        DO I=2,ILES
         AUX(I,K)=1./(AA(K)+CC(K)-BT(I-1,J-1)-CC(K)*EE(I,K-1))
         EE(I,K)=AA(K)*AUX(I,K)
         FF(I,K)=(-F(I,J,K)+CC(K)*FF(I,K-1))*AUX(I,K)
        ENDDO
  320  CONTINUE
       IF (J .EQ. 2) THEN
         Y1(2)=0.
       ELSE
         Y1(2)=(-F(2,J,KLES)+CC(KLES)*FF(2,KL2))
     1         /((1.-EE(2,KL2))*CC(KLES)-BT(1,J-1))
       ENDIF
       DO I=3,ILES
         Y1(I)=(-F(I,J,KLES)+CC(KLES)*FF(I,KL2))
     1         /((1.-EE(I,KL2))*CC(KLES)-BT(I-1,J-1))
       ENDDO
       DO I=2,ILES
        F(I,J,KLES)=Y1(I)
       ENDDO
       DO 390 K=3,KLES
         K1=KMAX+1-K
         DO I=2,ILES
           Y1(I)=EE(I,K1)*Y1(I)+FF(I,K1)
           F(I,J,K1)=Y1(I)
         ENDDO
  390  CONTINUE
 3000 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ***    BACKWARD TRANSFORMATION IN X AND Y    *********************
C     ***    BACKWARD TRANSFORMATION IN X

c$doacross local(I,IM,J,JM,K)
        DO I=2,IHF1
        IM=I-1
        DO J=2,JLES
        JM=J-1
        DO K=2,KLES
         AR(K-1,JM,IM)=F(I,J,K)
        ENDDO
        ENDDO
        ENDDO

c$doacross local(I,IM,J,JM,K)
        DO I=IHF2,ILES
        IM=I-1
        DO J=2,JLES
        JM=J-1
        DO K=2,KLES
         AI(K-1,JM,IM)=-F(I,J,K)
        ENDDO
        ENDDO
        ENDDO
      DO J=1,JL2
      DO K=1,KL2
       AI(K,J,1)=0.0
       AI(K,J,IHF)=0.0
      ENDDO
      ENDDO

c$doacross local(I,IR,J,K)
      DO 53 I=IHF1,IL2
        IR=IMAX-I
        DO J=1,JL2
        DO K=1,KL2
         AR(K,J,I)=AR(K,J,IR)
         AI(K,J,IR)=-AI(K,J,I)
        ENDDO
        ENDDO
   53 CONTINUE
C     ******      ******
      CALL FFTX (1,IL2)
C     ******      ******
c$doacross local(K,KM,J,JM,I)
      DO K=2,KLES
         KM=K-1
         DO J=2,JLES
            JM=J-1
            DO I=2,ILES
               F(I,J,K)=AR(KM,JM,I-1)
            ENDDO
         ENDDO
      ENDDO
C     ***    BACKWARD TRANSFORMATION IN Y  ****************************
      IF (IL2 .EQ. JL2) THEN

c$doacross local(J,JM,I,IM,K)
         DO J=2,JHF1
         JM=J-1
         DO I=2,ILES
         IM=I-1
         DO K=2,KLES
           AR(K-1,IM,JM)=F(I,J,K)
         ENDDO
         ENDDO
         ENDDO

c$doacross local(j,jm,i,im,k)
         DO J=JHF2,JLES
         JM=J-1
         DO I=2,ILES
         IM=I-1
         DO K=2,KLES
           AI(K-1,IM,JM)=-F(I,J,K)
         ENDDO
         ENDDO
         ENDDO
        DO I=1,IL2
        DO K=1,KL2
          AI(K,I,1)=0.0
          AI(K,I,JHF)=0.0
        ENDDO
        ENDDO

c$doacross local(j,jr,i,k)
       DO 63 J=JHF1,JL2
          JR=JMAX-J
         DO I=1,IL2
         DO K=1,KL2
           AR(K,I,J)=AR(K,I,JR)
           AI(K,I,JR)=-AI(K,I,J)
         ENDDO
         ENDDO
   63  CONTINUE
C     ******      ******
       CALL FFTY (1,JL2)
C     ******      ******
c$doacross local(k,km,j,jm,i)
       DO K=2,KLES
       KM=K-1
       DO J=2,JLES
       JM=J-1
       DO I=2,ILES
         F(I,J,K)=AR(KM,I-1,JM)
       ENDDO
       ENDDO
       ENDDO
      ELSE
C     ***    BACKWARD TRANSFORMATION IN Y  *****************************
        CALL SLVBY (JHF,JHF1,JHF2)
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
      SUBROUTINE FFTX (IFLG,L)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****   FFT FOR X-DIREXTION                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NXI=NX-2,NYJ=NY-2,NZK=NZ-2)
      PARAMETER (NIQ=NXI/4+1,NZKNYJ=NZK*NYJ)
      PARAMETER (ND1=NX*NY*NZ-NXI*NYJ*NZK-3*NYJ*NZK)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BSAT/ AR(NZKNYJ,NXI),Y1(NZKNYJ),Y2(NZKNYJ),Y5(NZKNYJ),
     1    y7(nd1)
      COMMON/BSAT1/ AI(NZKNYJ,NXI),Y3(NZKNYJ),Y4(NZKNYJ),Y6(NZKNYJ),
     1    y8(nd1)
      COMMON/RFFTX/ S(NIQ),IV(NXI)
      save
      JKMAX=KL2*JL2
      IF (IFLG .EQ. 1) GO TO 9
C   ******      ********************************************************
      DO I=1,NIQ
       S(I)=0.
      ENDDO
      DO I=1,NXI
       IV(I)=0
      ENDDO
      IE=0
      I=1
    6 IE=IE+1
      I=I+I
      IF(I.LT.L) GO TO 6
       PINT=8.*ATAN(1.)/FLOAT(L)
       LIMIT=L/4+1
       DO I=1,LIMIT
         S(I)=SIN((I-1)*PINT)
       ENDDO
      IV(1)=0
      DO 5 I=2,L
        IVI=0
        I1=I-1
        I2=I1/2
        DO 8 J=1,IE
         IVI=IVI+IVI+I1-I2-I2
         I1=I2
         I2=I2/2
    8   CONTINUE
        IV(I)=IVI
    5 CONTINUE
      RETURN
    9 LO2=L/2
      LO4=L/4
      N2=LO2
      DO 100 I=1,N2
        IR=I+N2
       DO KJ=1,JKMAX
         Y1(KJ)=AR(KJ,I)+AR(KJ,IR)
         Y2(KJ)=AI(KJ,I)+AI(KJ,IR)
         Y3(KJ)=AR(KJ,I)-AR(KJ,IR)
         Y4(KJ)=AI(KJ,I)-AI(KJ,IR)
         AR(KJ,I)=Y1(KJ)
         AI(KJ,I)=Y2(KJ)
         AR(KJ,IR)=Y3(KJ)
         AI(KJ,IR)=Y4(KJ)
       ENDDO
  100 CONTINUE
      MM=L
      JJ=1
   50 MM=MM/2
      N2=N2/2
      JJ=JJ+JJ
      KS=-1

      DO J=1,JJ
         KS=KS+1
         IX=IV(KS+KS+1)
         JMM=MM*(J-1)
         IF(LO4.GE.IX)THEN 
            SINV=S(IX+1)
            COSV=S(LO4-IX+1)
         ELSE
            SINV=S(LO2-IX+1)
            COSV=-S(IX-LO4+1)
         ENDIF
         DO I=1,N2
            K=I+JMM
            L2=K+N2
            DO KJ=1,JKMAX
               Y11=AR(KJ,L2)*COSV-AI(KJ,L2)*SINV
               Y22=AR(KJ,L2)*SINV+AI(KJ,L2)*COSV
               Y33=AR(KJ,K)+Y11
               Y44=AI(KJ,K)+Y22
               Y55=AR(KJ,K)
               Y66=AI(KJ,K)
               AR(KJ,L2)=y55-Y11
               AI(KJ,L2)=y66-Y22
               AR(KJ,K)=Y33
               AI(KJ,K)=Y44
            ENDDO
         ENDDO
      ENDDO

      IF(N2.GE.2) GO TO 50
      DO 60 I=1,L
        K=IV(I)+1
       IF(I.GE.K) GO TO 60
         DO KJ=1,JKMAX
           Y1(KJ)=AR(KJ,K)
           Y2(KJ)=AI(KJ,K)
           Y3(KJ)=AR(KJ,I)
           Y4(KJ)=AI(KJ,I)
          AR(KJ,K)=Y3(KJ)
          AI(KJ,K)=Y4(KJ)
          AR(KJ,I)=Y1(KJ)
          AI(KJ,I)=Y2(KJ)
         ENDDO
   60 CONTINUE
      RETURN
      END

      SUBROUTINE FFTY (IFLG,L)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ****   FFT FOR Y-DIREXTION                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NXI=NX-2,NYJ=NY-2,NZK=NZ-2)
      PARAMETER (ND1=NX*NY*NZ-NXI*NYJ*NZK-2*NXI*NZK)
      PARAMETER (NJQ=NYJ/4+1,NZKNXI=NZK*NXI)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BSAT/ AR(NZKNXI,NYJ),Y1(NZKNXI),Y2(NZKNXI),Y5(ND1)
      COMMON/BSAT1/ AI(NZKNXI,NYJ),Y3(NZKNXI),Y4(NZKNXI),Y6(ND1)
      COMMON/RFFTY/ S(NJQ),IV(NYJ)
      save
      IKMAX=KL2*IL2
      IF(IFLG.EQ.1) GO TO 9
C   ******      ********************************************************
      DO J=1,NJQ
       S(J)=0.
      ENDDO
      DO J=1,NYJ
       IV(J)=0
      ENDDO
      IE=0
      I=1
    6 IE=IE+1
      I=I+I
      IF(I.LT.L) GO TO 6
       PINT=8.*ATAN(1.)/FLOAT(L)
       LIMIT=L/4+1
       DO I=1,LIMIT
          S(I)=SIN((I-1)*PINT)
       ENDDO
      IV(1)=0
      DO 5 I=2,L
        IVI=0
        I1=I-1
        I2=I1/2
        DO 8 J=1,IE
         IVI=IVI+IVI+I1-I2-I2
         I1=I2
         I2=I2/2
    8   CONTINUE
      IV(I)=IVI
    5 CONTINUE
      RETURN
    9 LO2=L/2
      LO4=L/4
      N2=LO2

      DO I=1,N2
         IR=I+N2
         DO KI=1,IKMAX
            Y11=AR(KI,I)+AR(KI,IR)
            Y22=AI(KI,I)+AI(KI,IR)
            Y33=AR(KI,I)-AR(KI,IR)
            Y44=AI(KI,I)-AI(KI,IR)
            AR(KI,I)=Y11
            AI(KI,I)=Y22
            AR(KI,IR)=y33
            AI(KI,IR)=y44
         ENDDO
      ENDDO

      MM=L
      JJ=1
   50 MM=MM/2
      N2=N2/2
      JJ=JJ+JJ
      KS=-1

      DO 300 J=1,JJ
         KS=KS+1
         IX=IV(KS+KS+1)
         JMM=MM*(J-1)
         IF(LO4.GE.IX) THEN
            SINV=S(IX+1)
            COSV=S(LO4-IX+1)
         ELSE
            SINV=S(LO2-IX+1)
            COSV=-S(IX-LO4+1)
         ENDIF
         DO 30 I=1,N2
            K=I+JMM
            L2=K+N2
            DO KI=1,IKMAX
               Y11=AR(KI,L2)*COSV-AI(KI,L2)*SINV
               Y22=AR(KI,L2)*SINV+AI(KI,L2)*COSV
               Y33=AR(KI,K)+Y11
               Y44=AI(KI,K)+Y22
               AR(KI,L2)=AR(KI,K)-Y11
               AI(KI,L2)=AI(KI,K)-Y22
               AR(KI,K)=Y33
               AI(KI,K)=Y44
            ENDDO
   30    CONTINUE
  300 CONTINUE

      IF(N2.GE.2) GO TO 50
       DO 60 I=1,L
         K=IV(I)+1
        IF(I.GE.K) GO TO 60
         DO KI=1,IKMAX
           Y1(KI)=AR(KI,K)
           Y2(KI)=AI(KI,K)
          AR(KI,K)=AR(KI,I)
          AI(KI,K)=AI(KI,I)
          AR(KI,I)=Y1(KI)
          AI(KI,I)=Y2(KI)
        ENDDO
   60 CONTINUE
      RETURN
      END
      SUBROUTINE SLVFY (JHF1,JHF2,RJL2)
C     ****   FORWARD TRANSFORMATION IN Y   *****************************
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NXI=NX-2,NYJ=NY-2,NZK=NZ-2)
      PARAMETER (ND1=NX*NY*NZ-NXI*NYJ*NZK)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BSAT/ AR(NZK,NXI,NYJ),XXD(ND1)
      COMMON/BSAT1/ AI(NZK,NXI,NYJ),YYD(ND1)
      COMMON/BADV/ F(NX,NY,NZ)
      save
c$doacross local(k,km,j,jm,i)
      DO K=2,KLES
      KM=K-1
      DO J=2,JLES
      JM=J-1
      DO I=2,ILES
        AR(KM,I-1,JM)=F(I,J,K)*RJL2
      ENDDO
      ENDDO
      ENDDO

c$doacross local(j,i,k)
      do 25 j=1,jl2
        do i=1,il2
        do k=1,kl2
          ai(k,i,j)=0.0
        enddo
        enddo
   25 CONTINUE
C     ******      ******
      CALL FFTY (1,JL2)
C     ******      ******

c$doacross local(k,km,j,jm,i)
      DO K=2,KLES
        KM=K-1
        DO J=2,JHF1
          JM=J-1
          DO I=2,ILES
            F(I,J,K)=AR(KM,I-1,JM)
          ENDDO
        ENDDO
      ENDDO

c$doacross local(k,km,j,jm,i)
      DO K=2,KLES
        KM=K-1
        DO J=JHF2,JLES
          JM=J-1
          DO I=2,ILES
            F(I,J,K)=AI(KM,I-1,JM)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE SLVBY (JHF,JHF1,JHF2)
C     ****   BACKWARD TRANSFORMATION IN Y   ****************************
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NXI=NX-2,NYJ=NY-2,NZK=NZ-2)
      PARAMETER (ND1=NX*NY*NZ-NXI*NYJ*NZK)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BSAT/ AR(NZK,NXI,NYJ),XXD(ND1)
      COMMON/BSAT1/ AI(NZK,NXI,NYJ),YYD(ND1)
      COMMON/BADV/ F(NX,NY,NZ)
      save

c$doacross local(j,jm,i,im,k,km)
      DO J=2,JHF1
         JM=J-1
         DO I=2,ILES
            IM=I-1
            DO K=2,KLES
               KM=K-1
               AR(KM,IM,JM)=F(I,J,K)
            ENDDO
         ENDDO
      ENDDO

c$doacross local(j,jm,i,im,k,km)
      DO J=JHF2,JLES
         JM=J-1
         DO I=2,ILES
            IM=I-1
            DO K=2,KLES
               KM=K-1
               AI(KM,IM,JM)=-F(I,J,K)
            ENDDO
         ENDDO
      ENDDO

       do k=1,kl2
       do i=1,il2
         AI(k,i,1)=0.0
         AI(k,i,JHF)=0.0
       enddo
       enddo

      DO 63 J=JHF1,JL2
       JR=JMAX-J
       do i=1,nxi
       do k=1,nzk
         AR(k,i,J)=AR(k,i,JR)
         AI(k,i,JR)=-AI(k,i,J)
       enddo
       enddo
   63 continue
C     ******      ******
      CALL FFTY (1,JL2)
C     ******      ******
c$doacross local(k,km,j,jm,i)
      DO K=2,KLES
      KM=K-1
      DO J=2,JLES
      JM=J-1
      DO I=2,ILES
       F(I,J,K)=AR(KM,I-1,JM)
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
      subroutine bndop (x,x1)
c     ***  set cyclic lateral boundary condition   *********************
      parameter (NX=130,NY=130)
      common/bxyz/ imax,iles,il2,jmax,jles,jl2,kmax,kles,kl2,nisec(5)
      dimension x(nx,ny,1),x1(nx,ny,1)
      save
c     ***  set periodic boundary condition in x-axis
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 21 k=2,kles
      do 20 j=2,jles
       x(1,j,k)=x(iles,j,k)
       x1(1,j,k)=x1(iles,j,k)
       x(imax,j,k)=x(2,j,k)
       x1(imax,j,k)=x1(2,j,k)
   20 continue
   21 continue
c     ***  periodic boundary condition in y-axis
      do 41 k=2,kles
      do 40 i=1,imax
       x(i,1,k)=x(i,jles,k)
       x1(i,1,k)=x1(i,jles,k)
       x(i,jmax,k)=x(i,2,k)
       x1(i,jmax,k)=x1(i,2,k)
   40 continue
   41 continue
c     ***  for lower and upper boundary condition
c$doacross local(j,i)
      do 61 j=1,jmax
      do 60 i=1,imax
       x(i,j,1)=x(i,j,2)
       x1(i,j,1)=x1(i,j,2)
       x(i,j,kmax)=x(i,j,kles)
       x1(i,j,kmax)=x1(i,j,kles)
   60 continue
   61 continue
      return
      end

      SUBROUTINE SEPCA
CC   ***  CHURCHILL AND HOUZES' CONVECTIVE AND ANVIL SEPARATION   *****
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880)
      PARAMETER (NZ2=2*NZ,NZ7=7*NZ,NZ15=15*NZ+NX)
      PARAMETER (NT2=2*NT,NM13=13*NM)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5(4),AL,CP,RA,CK,CE,EPS(6)
      COMMON/B2TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B2CR/ QC(NX,NY,NZ),QR(NX,NY,NZ)
      COMMON/B2IG/ QI(NX,NY,NZ),QG(NX,NY,NZ)
      COMMON/B2S/ QS(NX,NY,NZ)
      COMMON/B4WP/ W(NX,NY,NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ15)
      COMMON/B6/ FD(NZ2),P0(NZ),PI(NZ),F0(NZ7)

      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NX,NY),ICS5(NX,NY,4),
     1  IBZ(NX,NY,4)
      COMMON/BI/ IT(NX,NY),ICS(NX,NY,4)
      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      COMMON/RSTAT/ CSTTT(NX,NY),CSTT(NX,NY)

      COMMON/BBA/ XYP(NX,NY),XY(NX,NY),XYM(NX,NY),XYPM(NX,NY),Y1(NM),
     1   TAIR(NM),Y2(NM),Y3(NM13)
      DIMENSION JV(NX,NY)
      save
CCC
         JOUT=JMAX+1
         iout=imax+1
         C2D25=2./25.
c$doacross local(j,i)
        DO J=1,JMAX
        DO I=1,IMAX
         CSTT(I,J)=0.
         JV(I,J)=0
         IT(I,J)=0
        ENDDO
        ENDDO
c$doacross local(ist,j,i)
        DO IST=2,3
        DO J=1,JMAX
        DO I=1,IMAX
           ICS(I,J,IST)=0
        ENDDO
        ENDDO
        ENDDO
      DO 101 J=2,JLES
       JPK=J
        JP1=JPK+1
        JP2=JPK+2
        JM1=JPK-1
        JM2=JPK-2
         IF(JP1.EQ.JMAX) JP1=2
         IF(JP2.EQ.JMAX) JP2=2
         IF(JP2.EQ.JOUT) JP2=3
         IF(JM1.EQ.1) JM1=JLES
         IF(JM2.EQ.1) JM2=JLES
         IF(JM2.EQ.0) JM2=JL2
      DO 100 I=2,iles
       IF (RI(I,J) .LT. 0.001) GOTO 100
       IPK=I
        IP1=IPK+1
        IP2=IPK+2
        IM1=IPK-1
        IM2=IPK-2
         IF(iP1.EQ.iMAX) iP1=2
         IF(iP2.EQ.iMAX) iP2=2
         IF(iP2.EQ.iOUT) iP2=3
         IF(iM1.EQ.1) iM1=iLES
         IF(iM2.EQ.1) iM2=iLES
         IF(iM2.EQ.0) iM2=iL2
      RV=RI(IM2,JM2)+RI(IM1,JM2)+RI(IPK,JM2)+RI(IP1,JM2)+RI(IP2,JM2)
       RV=RV+RI(IM2,JM1)+RI(IM1,JM1)+RI(IPK,JM1)+RI(IP1,JM1)+RI(IP2,JM1)
       RV=RV+RI(IM2,JPK)+RI(IM1,JPK)+RI(IPK,JPK)+RI(IP1,JPK)+RI(IP2,JPK)
       RV=RV+RI(IM2,JP1)+RI(IM1,JP1)+RI(IPK,JP1)+RI(IP1,JP1)+RI(IP2,JP1)
       RV=RV+RI(IM2,JP2)+RI(IM1,JP2)+RI(IPK,JP2)+RI(IP1,JP2)+RI(IP2,JP2)
      RV=C2D25*RV
       IF(RI(IPK,JPK).GE.RV)THEN
        IT(IPK,JPK)=3
         IF(RI(IM1,JPK) .GE. 0.01) IT(IM1,JPK)=3
         IF(RI(IP1,JPK) .GE. 0.01) IT(IP1,JPK)=3
         IF(RI(IM1,JP1) .GE. 0.01) IT(IM1,JP1)=3
         IF(RI(IP1,JP1) .GE. 0.01) IT(IP1,JP1)=3
         IF(RI(IPK,JP1) .GE. 0.01) IT(IPK,JP1)=3
         IF(RI(IM1,JM1) .GE. 0.01) IT(IM1,JM1)=3
         IF(RI(IPK,JM1) .GE. 0.01) IT(IPK,JM1)=3
         IF(RI(IP1,JM1) .GE. 0.01) IT(IP1,JM1)=3
      ELSE
        IF(IT(IPK,JPK).NE.3) IT(IPK,JPK)=2
      ENDIF
  100 CONTINUE
  101 CONTINUE
C  ***  FIND THE LOCATION AND VALUE OF MAX SFC RAINFALL GRID POINT
        BIG=0.0
       DO 151 J=2,JLES
       DO 150 I=2,iles
        A1=RI(I,J)
        IF (A1 .LT. BIG) GO TO 150
c        IBIG=I
        BIG=A1
  150  CONTINUE
  151  CONTINUE

cc        RIANMAX=MAX(20., MIN(25., BIG))
        RIANMAX=20.
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,iles
         IF (RI(I,J) .GE. RIANMAX) IT(I,J)=3
         IF (RI(I,J) .LT. 0.001) IT(I,J)=1
       ENDDO
       ENDDO
C   ***  TAO & SIMPSON (1989) AND TAO ET AL (1993)   **************

C     ****   SEARCH FOR THE LARGEST VALUE OF W AND CLOUD WATER
      BIGW=0.
      DO K=2,KLES
      DO J=2,JLES
      DO I=2,ILES
       XA=ABS(W(I,J,K))
       BIGW=MAX (BIGW,XA)
      ENDDO
      ENDDO
      ENDDO
cc       WBIG=MAX(500., MIN(500., 0.5*BIGW), 0.25*BIGW)
c      BIGQC=0.
c      DO K=2,KLES
c      DO J=2,JLES
c      DO I=2,ILES
c         XA=QC(I,J,K)
c         BIGQC=MAX (BIGQC,XA)
c      ENDDO
c      ENDDO
c      ENDDO
c       QCBIG= MAX(1.00E-3, 0.5*BIGQC)
c       QCBIG1=MAX(2.00E-3, 0.5*BIGQC)
       wbig=300.
       qcbig= 0.50e-3
       qcbig1=1.00e-3
      do k=2,kles
         tair(k)=pi(k)*tb(k)-273.16
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO 351 J=2,JLES
         DO 350 I=2,iles
          ICLOUDW=0
          ICLOUDL=0
          ICLOUD=0
          ICWWL=0
          DO 35 K=2,KLES
c            TAIR(K)=PI(K)*(TA(K)+DPT(I,J,K))-273.16
CC  ***   LOWER CLOUDY REGION
           IF (TAIR(K) .GE. 0.0) THEN
            IF (QC(I,J,K) .GE. QCBIG) ICLOUDW=1
            IF (QC(I,J,K) .GE. QCBIG1) ICLOUDL=1
            IF (W(I,J,K) .GE. WBIG) ICWWL=1
           ENDIF
CC  ***   MIDDLE AND UPPER CLOUDY REGION
c           IF (TAIR(K) .LT. 0.0) THEN
c            IF (W(I,J,K) .GE. WBIG) ICWWU=1
c           ENDIF
            RZM=RHO(K)*1.E6
           Y3(K)=RZM*(QC(I,J,K)+QI(I,J,K)+QR(I,J,K)+QS(I,J,K)+QG(I,J,K))
           IF (Y3(K) .GE. .01) ICLOUD=1
   35     CONTINUE
C   ***  MIDDLE-UPPER LEVEL W > 0.5*WBIG M/S OR LOW-LEVEL QC > 1 G/KG
          IF (IT(I,J).EQ.1) THEN
            IF (ICWWL .EQ. 1 .AND. ICLOUDW .EQ. 1) IT(I,J)=3
          ENDIF
           IF (IT(I,J).EQ.2) THEN
             IF (ICLOUDL .EQ. 1) IT(I,J)=3
           ENDIF
C          IF (IT(I,J).EQ.2) THEN
C            IF (I .GT. IBIG-5) THEN
C             IF (ICWWU .EQ.1) THEN
C               IT(I,J)=3
C               IT(I+1,J)=3
C               IT(I+2,J)=3
C               IT(I+3,J)=3
C               IT(I+4,J)=3
C             ENDIF
C            ENDIF
C          ENDIF
          IF (IT(I,J).EQ.1 .AND. ICLOUD.EQ.1) JV(I,J)=1
  350   CONTINUE
  351   CONTINUE
C   ***   ********************************************************
C   ***   IT(I)=3  CONVECTIVE REGION
C   ***   IT(I)=2  STRATIFORM REGION
C   ***   IT(I)=4  STRATIFORM REGION BUT NO SFC PRECIPITATION
C   ***   IT(I)=1  CLOUD FREE REGION

       A1=.0
       A2=.0
       A3=.0
       A4=.0
       LCONV5=0
       LANVL5=0
       LNSPT5=0
       DO 900 J=2,JLES
        DO 90 I=2,iles
         IF(IT(I,J) .EQ. 2) THEN
          ICS(I,J,3)=1
	  ICS5(I,J,3)=1
          LS=LS+1
          A2=A2+RI(I,J)
          A3=A3+QR(I,J,2)
         ENDIF
   90   CONTINUE
        DO 92 I=2,iles
        IF(IT(I,J) .EQ. 3) THEN
          ICS(I,J,2)=1
          ICS5(I,J,2)=1
          LC=LC+1
          A1=A1+RI(I,J)
          A4=A4+QR(I,J,2)
          LCONV=LCONV+1
        ENDIF
   92   CONTINUE
        DO 94 I=2,iles
        IF(IT(I,J) .EQ. 1 .AND. JV(I,J) .EQ. 1) THEN
          ICS(I,J,4)=1
          ICS5(I,J,4)=1
          IT(I,J)=4
          LNSPT=LNSPT+1
        ENDIF
        IVV(I,J)=IT(I,J)
   94   CONTINUE
  900   CONTINUE
c$doacross local(j,i)
        DO 951 J=2,JLES
        DO 950 I=2,iles
          ICS5(I,J,4)=ICS(I,J,4)
          ICS5(I,J,2)=ICS(I,J,2)
          ICS5(I,J,3)=ICS(I,J,3)
          CSTT(I,J)=IT(I,J)
  950   CONTINUE
  951   CONTINUE
       ACO5=ACO5+A1
       AAN5=AAN5+A2
       ACO15=ACO15+A4
       AAN15=AAN15+A3
      RETURN
      END
ccc
cc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cc
ccc
      SUBROUTINE PBLIN (u,v,dpt,dqv,ri)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ****   SET DOMAIN AND INITIAL CONDITION
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (ILES=NX-1,JLES=NY-1)
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
C
      COMMON/DINRAD/ P00(NZ),DZ0(NZ),TAIRSFC(NX,NY)
      common/dinradn/ qairsfc(nx,ny),pairsfc(nx,ny),thairsf(nx,ny)
      COMMON/BPBL/ UHT(NZ),WHT(NZ),TGBAT0
      COMMON/BGWY/ UBAR,VBAR
      common/dinrad1/ dz1half

C      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
C 
      COMMON /SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
      DIMENSION UBOT(NX),VBOT(NX),TBOT(NX),QBOT(NX),RAINR(NX)
      DIMENSION RSHRT(NX),RLONG(NX)
      dimension u(nx,ny),v(nx,ny),dpt(nx,ny),dqv(nx,ny),ri(nx,ny)
      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      IF (IFLAG .EQ. 1) GOTO 2000
c      RETURN
C
c 2000 CONTINUE
C
       r100=1./100.
       r1000=1./1000.
      DO J=2,JLES
        DO I=2,ILES
        UBOT(I)=(UBAR+.5*(U(I,J)+U(I+1,J)))*r100
        VBOT(I)=(VBAR+.5*(V(I,J)+V(I,J+1)))*r100 
        TBOT(I)=(TA1(2)+DPT(I,J))*pi(2)-273.16  ! potential T
        QBOT(I)=1000.*(QA1(2)+DQV(I,J))         ! vapor mixing ratio
        RSHRT(I)=140.
        RLONG(I)=400.                           ! for now (does not need)
         RAINR(I)=RI(I,J)  
C       PRINT *, UBOT(I),VBOT(I),TBOT(I),QBOT(I)
        ENDDO    
C       HT=20.
c
       HT=dz1half
       SSTIN=tairsfc(i,j)-273.16
       PSFCMB=pairsfc(i,j)*r1000 ! pressure at surfase (from microbar to mb)
c       ATIME=.2                  ! for now
C       print *, sstin,psfcmb 

c       CALL SFFLUX (ATIME,UBOT,VBOT,HT,SSTIN,TBOT,J,
       CALL SFFLUX (UBOT,VBOT,HT,SSTIN,TBOT,J,
     1              QBOT,RSHRT,RLONG,RAINR,PSFCMB)
C
      ENDDO
C
      write(29) suw
      write(29) svw
      write(29) swt
      write(29) swq

      RETURN 
      END
 
      SUBROUTINE SFFLUX (UBOT,VBOT,HT,SSTIN,TBOT,JY,
     1                   QBOT,RSHRT,RLONG,RAINR,PSFC)
      PARAMETER (NX=130,NY=130)
      COMMON/SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
      DIMENSION UBOT(1),VBOT(1),TBOT(1),QBOT(1),RAINR(1)  
      DIMENSION RSHRT(1),RLONG(1)
      DIMENSION USTARs(NX),TSTARs(NX),QSTARs(NX)

c COARE bulk flux version 2.0 (Toulouse release)
c Based on Liu et al. 1979 and Liu's original computer code.
c 1   First modified by David Rogers October 12, 1993
c 2   Adopt blending between Kansas and Free-convection forms by Fairall,
c     Rogers and Bradley October 20 1993
c        Modified by Fairall and Bradley June 22 1994 to include
c 3   cool skin effect(needs radiative fluxes)
c     equations quoted are from Fairall,Bradley and Godfrey (unpub ms 1994)
c     NB if IR SST used, set Jcool=0
c 4   sensible heat flux (ocean cooling) due to rain at wet bulb temp.
c     formalism by Gosnel,Webster and Fairall (unpub ms 1994)
c 5   a simplified version of Price, Weller and Pinkel (PWP) model for solar warm layer
c     this option requires a long data file in time order, 
c     with resolution at least 1 hour
c 6   Subroutine H_ADJUST added by Meghan Cronin 6/25/94 to adjust ws,qq,Tair
c     to specified heights (e.g. 10 m). Wind speed can be adjusted to standard
c     height different than the humidity and temperature standard height. 
c....................................................................
c   input:
c     Intime (time - fraction of a day) days          SPA       
c     hUm (wind measurement height) m       MC
c     hTm (T&rh measurement height) m       MC
c     hUs (wind standard height) m       MC
c     hTs (T&rh standard height) m       MC
c     ts_depthx (depth of sst instrument) m 
c     ws (wind speed) m/s
c     sst (sea surface temp.)  deg. C
c     atb (air temperature) deg. C
c     qq (specific humidity) g/kg or qq (RH as decimal) NB Code works in kg/kg!!
c     rs (shortwave radiation) W/m2
c     rl (downwelling longwave) W/m2
c     rain (average rainrate) mm/hour
c     pp (pressure) mb
c     zi (boundary-layer depth; use 600m if data unavailable) m
c     Jcool (=1 for cool skin calculation; =0 if SST by IR radiometer)
c     Jwarm (=0 to ignore; otherwise =2 in first line of data and =1 all other lines
c
c   output:
c
c     HF W/m**2   (Sensible heat flux)
c     EF W/m**2   (Latent heat flux)
c     RF W/m**2   (Rainfall heat flux)
c     TAU m**2/s**2
c     Ustar m/s
c     Qstar kg/kg
c     Tstar C
c     CD - drag coefficient
c     CH - transfer coefficient for heat
c     CE - transfer coefficient for moisture
c     RR - Roughness Reynolds number
c     RT - Roughness Reynolds number for temperature
c     RQ - Roughness Reynolds number for moisture
c     ZL - ht/L where L is the Obukhov length
c     Z0 - roughness length
c     zot - roughness length for temperature
c     zoq - roughness length for humidity
c     T0 - skin temperature C
c

      real*8 hUm,hTm,hUs,hTs,ws_h,qq_h,ta_h
c      real*8 loc,time
      real*8 ws,sst,atb,qq,pp,zi,rain
      real*8 QH,QE,TAU,Ustar,Qstar,Tstar
      real*8 rl,rs,RF,T0
c     real*8 ts_depth
      real*8 CD,CE,CH,RR,RT,RQ,Zl,ZO
      real*8 Jcool,Jwarm,zot,zoq,dt_wrm,dter
      real*8 jtime,qcol_ac,tau_ac
      integer jamset
      common /old/jtime,qcol_ac,tau_ac,tau_old,rf_old,hf_old,ef_old,
     1            jamset
      save
c
c
c      loc=10                     ! local offset from GMT
c ---- housekeeping variables --------
      qcol_ac=0.
      tau_ac=0.
      jtime=0
      jamset=0
      tau_old=0.
      hf_old=0.
      ef_old=0.
c--------------------------------------
C      
c ---------- read in data -------------
c      open(unit=10,file='input_data',status='old')
c      open(unit=12,file='output_data')
c if instruments are at fixed heights they can be read in or set here
      hUm=HT
      hTm=HT
      hUs=10.                  ! standard reference height
      hTs=10.
c      ts_depth=.05
c 
c set jcool = 1 for cool skin calculation or = 0 if sst from IR radiometer
      jcool=0
c 
c set jwarm = 2 for warm layer calc., = 0 to ignore or if sst from IR radiometer
      jwarm=0
c
c read in time series of data sequentially
c
c test data
c      time=ATIME   ! fraction of a day (<=1)
      zi=600.       ! m
      pp=PSFC       ! mb
      sst=SSTIN     ! degrees C

      do i=2,NX-1
       ws= SQRT(UBOT(I)**2+VBOT(I)**2)
       atb=TBOT(I)                          ! degrees C
       qq=QBOT(I)                           ! g/kg
       rs=RSHRT(I)                          ! W/m^2
       rl=RLONG(I)                          ! W/m^2
       rain=RAINR(I)                        ! mm/hr
c1      read(10,*,end=999) time,ws,sst,atb,qq,rs,rl,rain,pp,zi
       call bulk_flux (hUm,hTm,hUs,hTs,ws,sst,atb,
     &   qq,WS_H,TA_H,QQ_H,rs,rl,
     &   rain,pp,zi,jcool,jwarm,QH,QE,RF,TAU,USTAR,TSTAR,QSTAR,CD,CH,
c     &   CE,RR,RT,RQ,ZL,ZO,ZOT,ZOQ,DT_WRM,DTER,T0,ts_depth)
     &   CE,RR,RT,RQ,ZL,ZO,ZOT,ZOQ,DT_WRM,DTER,T0)
       if (i .eq. 600) then
        print *, 'sensible heat flux, latent heat flux, ustar'
        print *, qh, qe, ustar   
       endif
c
C      write(*,*) time,ws,sst,TA_H,qq_H,T0,rs,rl,rain,pp,
C    &  QH,QE,RF,TAU,USTAR,TSTAR,QSTAR,CD,CH,CE,ZL,ZO,DT_WRM,DTER
       ustarS(I)=ustar
       tstarS(I)=tstar
       qstarS(I)=qstar
      enddo
C
C
      do i=1,nx
       suw(i,JY)=0.
       svw(i,JY)=0.
       swt(i,JY)=0.
       swq(i,JY)=0.
      enddo
C
C        ws= SQRT(UBOT(2)**2+VBOT(2)**2)            
C        suw(1,JY)=-(2.*ustars(2)-ustars(3))**2*ubot(2)/ws 
C       svw(1,JY)=-(2.*ustars(2)-ustars(3))**2*vbot(2)/ws
C       SWT(1,JY)=-USTARS(2)*TSTARS(2)
C       SWQ(1,JY)=-USTARS(2)*QSTARS(2)


      do i=2,nx-1
       WS= SQRT(UBOT(I)**2+VBOT(I)**2)
       SUW(I,JY)=-USTARS(I)**2*UBOT(I)/WS
       SVW(I,JY)=-USTARS(I)**2*VBOT(I)/WS
       SWT(I,JY)=-USTARS(I)*TSTARS(I)
       SWQ(I,JY)=-USTARS(I)*QSTARS(I)
      enddo
      
C     write(29) suw
C     write(29) svw
C     write(29) swt
C     write(29) swq

      return
      end

      subroutine bulk_flux(
     & hUm,hTm,hUs,hTs,                       !    MC
     & ws,sst,atb,qq,ws_h,Ta_h,qq_h,          !    MC
     & rs,rl,rainx,pp,zix,Jcoolx,Jwarmx,HF,EF,RF,TAU,
     & Ustar,Tstar,Qstar,
     & CD,CH,CE,RRx,RTx,RQx,ZLx,ZOx,zotx,zoqx,
c     & dt_wrmx,dterx,T0,ts_depthx)
     & dt_wrmx,dterx,T0)
      real*8 hUm,hTm,hUs,hTs,ws_h,Ta_h,qq_h
c      real*8 ts_depthx
      real*8 ws,sst,atb,qq,pp,zix,rainx
      real*8 HF,EF,TAU,Ustar,Qstar,Tstar
      real*8 rl,rs,RF,T0
      real*8 CD,CE,CH,RRx,RTx,RQx,Zlx,ZOx
      real*8 Jcoolx,Jwarmx,zotx,zoqx,dt_wrmx,dterx
      real*8 jtime,qcol_ac,tau_ac
      common /old/jtime,qcol_ac,tau_ac,tau_old,rf_old,hf_old,ef_old
     1           ,jamset
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
      common/wgust/DU_Wg,Wg     !MC
      save
      Jcool=Jcoolx
      Jwarm=Jwarmx
c............. MC added
         ZU=hUm       !height of wind measurement
         ZT=hTm       !height of temperature measurement
         ZQ=hTm       !height of water vapor measurement
         ZUs=hUs      !standard height of wind measurement
         ZTs=hTs      !standard height of temperature measurement
         ZQs=hTs      !standard height of water vapor measurement
c............ MC end
         U=ws        !wind speed m/s
         TS=sst      !surface temp. Celsius
         T=atb       !air temp. Celsius
         P=pp        !pressure mb
         zi=zix
         toK=273.16
         Rnl= 0.97*(5.67e-8*(TS+toK)**4-rl)    ! Net longwave (up = +)
         Rns=0.945*rs                           ! Net shortwave (into water)
         rain=rainx
c        ts_depth=2.0           ! Franklin thermosalinograph intake
c         ts_depth=ts_depthx    ! depth of sst measurement
c -------------------- correct SST with PWP model ------------
      if(jwarm.eq.0) go to 15                   ! by-pass warm layer
Cwang    newtime=mod(loc+intime*24,24)*3600     !run time in secs SPA
Cwang if(jwarm.eq.2.) then                      ! first line of data
Cwang    jump=1                                 
Cwang    go to 16                               ! set jtime and pass thru' ASL
Cwang end if
Cwang if(newtime.gt.21600.and.jump.eq.1) goto 16 ! 6 am too late to start
Cwang if(newtime.lt.jtime) then                  ! reset at midnight
Cwang    jamset=0                                ! test threshold q morning only
Cwang    fxp=.5
Cwang    tk_pwp=19
Cwang    ts_pwp=ts
Cwang    tau_ac=0
Cwang    qcol_ac=0
Cwang    dt_wrm=0.
Cwang    jump=0
Cwang    rich=.65                                  ! critical Rich. No
Cwang    ctd1=sqrt(2*rich*cpw/(al*grav*rhow))        ! Chris integrates u*^2 so
Cwang    ctd2=sqrt(2*al*grav/(rich*rhow))/(cpw**1.5) ! has /rhoa in both of these
Cwang    go to 16
Cwang else
Cwang    dtime=newtime-jtime                       ! delta time
Cwang    qr_out=rnl+hf_old+ef_old+rf               ! flux out from previous pass
Cwang    q_pwp=fxp*rns-qr_out                      ! effective net warming
Cwang    if(q_pwp.lt.50.and.jamset.eq.0) goto 16 ! threshold to start integrating
Cwang    jamset=1
Cwang    tau_ac=tau_ac+tau_old*dtime          ! tau from previous pass
Cwang    if(qcol_ac+q_pwp*dtime.gt.0) then
Cwang      do index=1,5                       ! iterate for warm layer thickness
Cwang      fxp=1.-(0.28*0.014*(1-exp(-tk_pwp/0.014))
Cwang&       +0.27*0.357*(1-exp(-tk_pwp/0.357))
Cwang&       +.45*12.82*(1-exp(-tk_pwp/12.82)))/tk_pwp    ! solar absorb. prof
Cwang      qjoule=(fxp*rns-qr_out)*dtime
Cwang      if(qcol_ac+qjoule.gt.0.)             
Cwang&     tk_pwp=MIN(19.,ctd1*tau_ac/sqrt(qcol_ac+qjoule))
C          end do
C        else
C          fxp=.76
C          tk_pwp=19
C        endif
C       qcol_ac=qcol_ac+qjoule                  !integrate heat input
C       if(qcol_ac.gt.0) then
C         dt_wrm=ctd2*(qcol_ac)**1.5/tau_ac     ! pwp model warming
C       else
C         dt_wrm=0.
C       endif         
C     endif
C     if(tk_pwp.lt.ts_depth) then               ! pwp layer deeper than sensor
C       dsea=dt_wrm
C     else
C       dsea=dt_wrm*ts_depth/tk_pwp             ! linear temperature profile
C     endif
C       ts=ts+dsea        
C16    jtime=newtime
c--------------------------------- end warm layer ------------------------
15       call humidity(T,P,QA)      !Teten's formula returns sat. air in mb
      if(qq.lt.2.) then             !checks whether humidity in g/Kg or RH      
         R=qq
         ee=QA*R                    !convert from RH using vapour pressure      
         Q=.62197*(ee/(P-0.378*ee)) ! Spec. humidity kg/kg
      else
         Q=qq/1000.                 !g/kg to kg/kg
      endif
       QA=.62197*(QA/(P-0.378*QA)) !convert from mb to spec. humidity  kg/kg
       call humidity(TS,P,QS)        !sea QS returned in mb      
       QS=QS*0.98                    !reduced for salinity Kraus 1972 p. 46
       QS=.62197*(QS/(P-0.378*QS)) !convert from mb to spec. humidity  kg/kg
c
c----------------------calculate atmospheric surface layer ----------------
         CALL ASL(Jcool,IER)
         IF(IER.ge.0) then
C
C     COMPUTE SURFACE STRESS TAU, SENSIBLE HEAT FLUX HF,  
C     LATENT HEAT FLUX EF & other parameters
C
            TAU=rhoa*USR*usr*u/sqrt(u*u+wg*wg)
            HF=-cpa*rhoa*USR*TSR
            EF=-hl*rhoa*USR*QSR
c
c     compute heat flux due to rainfall
c
         dwat=2.11e-5*((T+toK)/toK)**1.94            ! water vapour diffusivity
         dtmp=(1.+3.309e-3*T-1.44e-6*T*T)*0.02411/(rhoa*cpa)! heat diffusivity
         dqs_dt=QA*hl/(rgas*(T+toK)**2)                     ! Clausius-Clapeyron
         alfac= 1/(1+0.622*(dqs_dt*hl*dwat)/(cpa*dtmp))     ! wet bulb factor
         RF= rain*alfac*cpw*((TS-T)+(QS-Q)*hl/cpa)/3600.
c
c -------------------------------- cool skin parameters ---------------
c
        dterx=dter
        T0=ts-dter
        tau_old=tau
        ef_old=ef
        hf_old=hf
c-------------------------------- warm layer parameter ----------------
CWANG   dt_wrmx=dt_wrm  
        dt_wrmx=0.      
C
C       COMPUTE TRANSFER COEFFICIENT
C

          CD=(USR/U)**2
          CH=USR*TSR/(U*(T-TS+.0098*zt+dter))
          CE=USR*QSR/(U*(Q-QS+dqer))                                      
            Ustar=USR
            Tstar=TSR
            Qstar=QSR
            RRx=RR
            RTx=RT
            RQx=RQ
            ZLx=ZL
            ZOx=ZO
            zotx=zot
            zoqx=zoq
c........... MC added
            ihumid=2
            if(qq .lt. 2) ihumid=1
               call h_adjust(ZUs,ZTs,ZQs,U_hs,T_hs,Q_hs,ihumid)
            ws_h=U_hs
            Ta_h=T_hs
            qq_h=Q_hs    
C......... MC end
         else
c input parameters out of range
            EF=-999.
            HF=-999.
            TAU=-999.
            Ustar=-999.
            Tstar=-999.
            Qstar=-999.
            RRx=-999.
            RTx=-999.
            RQx=-999.
            ZLx=-999.
            ZOx=-999.
c......... MC added
            ws_h=-999.
            Ta_h=-999.
            qq_h=-999.
c........... MC end
         endif
      return
      end
c ---------------------------------------------------------------------
      SUBROUTINE ASL(Jcool,IER)
c
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
      common/wgust/DU_Wg,Wg      ! MC
      save
c
C       TO EVALUATE SURFACE FLUXES, SURFACE ROUGHNESS AND STABILITY OF
C       THE ATMOSPHERIC SURFACE LAYER FROM BULK PARAMETERS BASED ON
C       LIU ET AL. (79) JAS 36 1722-1735 
C
c---------------------------  Factors  -------------------------------
         Beta=1.2     ! evaluated from Fairalls low windspeed turbulence data
         Von=0.4      ! von Karman's "constant"
         fdg=1.00     ! Fairall's LKB rr to von karman adjustment
         toK=273.16   ! Celsius to Kelvin
         grav=9.72    ! gravity equatorial value (ref. IGPP-SIO)
c--------------------------- Air constants ---------------------------
         Rgas=287.1                  ! J/kg/K     gas const. dry air
         hl=(2.501-0.00237*TS)*1e+6  ! J/kg  latent heat of vaporization at TS
         Cpa=1004.67           ! J/kg/K specific heat of dry air (Businger 1982)
c         Cpv=Cpa*(1+0.84*Q)     ! Moist air - currently not used (Businger 1982)
         rhoa=P*100./(Rgas*(T+toK)*(1.+.61*Q)) ! kg/m3  Moist air density ( "  )
         visa=1.326e-5*(1+6.542e-3*T+8.301e-6*T*T-4.84e-9*T*T*T)   ! m2/s
c      Kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
c--------------------------- Cool skin constants ---------------------
         al=2.1e-5*(ts+3.2)**0.79     ! water thermal expansion coefft.
         be=0.026                     ! salinity expansion coefft.
         cpw=4000.                    ! J/kg/K specific heat water
         rhow=1022.                  ! kg/m3  density water
         visw=1.e-6                   ! m2/s kinematic viscosity water
         tcw=0.6                      ! W/m/K   Thermal conductivity water
       bigc=16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)
       wetc=.622*hl*QS/(rgas*(TS+toK)**2) ! correction for dq;slope of sat. vap.
c
c---------------------------- Initialise everything  ---------------------
       IER=0
       ZL=0.                         
       Dter=0.                        ! cool skin Dt
       Dqer=0.                        ! cool skin Dq
c----------------------------  Initial guesses  ---------------------------
       US=0.                         !surface current = 0.
       Wg=0.5                        !Gustiness factor initial guess
       ZO=.0005                      !roughness initial guess
       tkt=.001                      ! guess sublayer thickness
       DU=U-US
       DU_Wg=(DU**2.+Wg**2.)**.5     !include gustiness in wind spd. difference
       DT=T-TS+.0098*zt              !potential temperature diff        
       DQ=Q-QS
       USR=.04*DU_Wg                 !
       TSR=.04*DT                    !initial guesses                        
       QSR=.04*DQ                    !
       IF(DU_Wg.NE.0.)THEN
       TA=T+toK
       RI=grav*ZU*(DT+0.61*TA*DQ)/(TA*DU_Wg**2)
       ELSE
       IER=-1
       RI=-999.
       ENDIF
       IF(RI.GT.0.25)IER=-1
c -----------------------------  Iterate 20 times  ------------------------
      do index=1,20
       CALL ZETA(T,Q,USR,TSR,QSR,ZU,ZLN)
       ZL=ZLN
       PUZ=PSI(1,ZL)
       ZTL=ZL*ZT/ZU
       ZQL=ZL*ZQ/ZU
       PTZ=PSI(2,ZTL)
       PQZ=PSI(2,ZQL)
       ZO=0.011*USR*USR/grav + 0.11*visa/USR        !after Smith 1988 
       USR=DU_Wg*von/(LOG(ZU/ZO)-PUZ)              !Gustiness effect incl.
       RR=ZO*USR/VISA
       CALL LKB(RR,RT,1)
       IF(RT.NE.-999.) GOTO 21
       IER = -2                                     !error - return
       RETURN
   21 CALL LKB(RR,RQ,2)
      IF(RQ.NE.-999.) GOTO 22
      IER = -2                                      !error - return
      RETURN
   22 zot=rt*visa/usr
      zoq=rq*visa/usr
      S = (LOG(ZT/zot)-PTZ)/(von*fdg)       !coeff fdg=1.04 included following
      D = (LOG(ZQ/zoq)-PQZ)/(von*fdg)       !Fairall observations during COARE. 
                                             !NOTE coeff changed to 1.
      dtt=(dt+dter)
      dqq=(dq+dqer)                         !or we could calculate new sat. hum.
      tsr=dtt/S                              !! modify
      qsr=dqq/D                              !! fluxes
      TVSR=TSR*(1.+0.61*Q)+(0.61*TA*QSR)
      Bf=-grav/TA*USR*TVSR
      if(Bf.gt.0) then
        Wg=Beta*(Bf*zi)**0.333  
      else
        Wg=0.
      endif
      DU_Wg=(DU**2.+Wg**2.)**.5                  !include gustiness in wind spd.
c --------------------------------  Cool skin  ---------------------------------
      if(Jcool.eq.0) goto 24
      hsb=-rhoa*cpa*usr*tsr
      hlb=-rhoa*hl*usr*qsr
      qout=rnl+hsb+hlb
      dels=rns*(.137+11*tkt-6.6e-5/tkt*(1-exp(tkt/8.0e-4))) ! Eq.16 Shortwave
      qcol=qout-dels
      if(qcol.gt.0.) then
        alq=Al*qcol+be*hlb*cpw/hl                     ! Eq. 7 Buoy flux water
        xlamx=6/(1+(bigc*alq/usr**4)**.75)**.333      ! Eq 13 Saunders coeff.
        tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)          ! Eq.11 Sublayer thickness
        dter=qcol*tkt/tcw                             ! Eq.12 Cool skin
      else
        dter=0.
      endif    
        dqer=wetc*dter
   24 continue
c--------------------------------- End cool skin  -------------------------------
      end do
c---------------------------------  End iterations  -----------------------------
      RETURN                                         ! to bulk_flux
      END
c
c
      Subroutine humidity(T,P,Qsat)                                 
      save
c
c     Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 
c
      Qsat = (1.0007+3.46e-6*P)*6.1121*exp(17.502*T/(240.97+T))
      return
      end
c
c
      SUBROUTINE LKB(RR,RT,IFLAG)
C
C       TO DETERMINE THE LOWER BOUNDARY VALUE RT OF THE 
C       LOGARITHMIC PROFILES OF TEMPERATURE (IFLAG=1) 
C       OR HUMIDITY (IFLAG=2) IN THE ATMOSPHERE FROM ROUGHNESS 
C       REYNOLD NUMBER RR BETWEEN 0 AND 1000.  OUT OF RANGE
C       RR INDICATED BY RT=-999. BASED ON LIU ET AL.(1979)
C       JAS 36 1722-1723
c     New scalar RR relation from Moana Wave data.
C
      DIMENSION A(9,2),B(9,2),RAN(9)
      save
      DATA A/0.177,2.7e3,1.03,1.026,1.625,4.661,34.904,1667.19,5.88E5,
     10.292,3.7e3,1.4,1.393,1.956,4.994,30.709,1448.68,2.98E5/
      DATA B/0.,4.28,0,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
     10.,4.28,0,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
      DATA RAN/0.11,.16,1.00,3.0,10.0,30.0,100.,300.,1000./
      I=1
      IF(RR.LE.0..OR.RR.GE.1000.) GOTO 90
   10 CONTINUE
      IF(RR.LE.RAN(I)) GOTO 20
      I=I+1
      GOTO 10
   20 RT=A(I,IFLAG)*RR**B(I,IFLAG)
      GOTO 99
   90 RT=-999.
   99 RETURN
      END
c
c
      FUNCTION PSI(ID,ZL)
C
C       TO EVALUATE THE STABILITY FUNCTION PSI FOR WIND SPEED (IFLAG=1)
C       OR FOR TEMPERATURE AND HUMIDITY PROFILES FROM STABILITY 
C       PARAMETER ZL. SEE LIU ET AL (1979).
c       Modified to include convective form following Fairall (Unpublished)
C
c      IF(ZL)10,20,30                                                 
       IF(ZL.lt.0.)goto 10                                                 
       IF(ZL.gt.0.)goto 30                                                 
       IF(ZL.eq.0.)goto 20                                                 
   10   F=1./(1+zl*zl) 
        CHIK=(1.-16.*ZL)**0.25
         IF(ID.EQ.1) GOTO 11
        PSIK=2.*LOG((1.+CHIK*CHIK)/2.)
        GOTO 12
   11   PSIK=2.*LOG((1.+CHIK)/2.)+LOG((1.+CHIK*CHIK)/2.)
     1 -2.*ATAN(CHIK)+2.*ATAN(1.)
   12   CHIC=(1.-12.87*ZL)**.333    !for very unstable conditions
        PSIC=1.5*LOG((CHIC*CHIC+CHIC+1.)/3.)
     &      -(3.**.5)*ATAN((2*CHIC+1.)/(3.**.5))
     &      +4.*ATAN(1.)/(3.**0.5) 
c                     
c     match Kansas and free-conv. forms with weighting F
c
      PSI= F*PSIK+(1-F)*PSIC
        goto 99
   20 PSI=0.
      GOTO 99
   30 continue
       PSI=-4.7*ZL
   99 RETURN
      END
      SUBROUTINE ZETA(T,Q,USR,TSR,QSR,Z,ZL)
C
C       TO EVALUATE OBUKHOVS STABILITY PARAMETER Z/L FROM AVERAGE
C       TEMP T IN DEG C, AVERAGE HUMIDITY Q IN GM/GM, HEIGHT IN M,
C       AND FRICTIONAL VEL,TEMP.,HUM. IN MKS UNITS
C       SEE LIU ET AL. (1979)
C
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
      save
c
      TA=T+toK
      TV=TA*(1.+0.61*Q)
      TVSR=TSR*(1.+0.61*Q)+0.61*TA*QSR
      IF(TVSR.EQ.0.)GOTO 10
      OB=TV*USR*USR/(grav*VON*TVSR)
      ZL=Z/OB
      GOTO 99
   10 ZL=0.
   99 RETURN
      END
       
c-------------------------------------------------------------------------
c
      SUBROUTINE H_ADJUST(ZUs,ZTs,ZQs,U_hs,T_hs,Q_hs,IHUMID)
C        This subroutine adjusts the U,T,Q variables to the specified
C        standard height (ZUs,ZTs,ZQs) using the loglayer profiles.
C        The DELTA correction (adjustment) is relative to the surface
C        measurement.             Cronin 4/13/94
C
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
      common/wgust/DU_Wg,Wg
      save
      CALL ZETA(T,Q,USR,TSR,QSR,ZUs,ZUsL)
      CALL ZETA(T,Q,USR,TSR,QSR,ZTs,ZTsL)
      CALL ZETA(T,Q,USR,TSR,QSR,ZQs,ZQsL)
      PUZs= PSI(1,ZUsL)
      PTZs= PSI(2,ZTsL)
      PQZs= PSI(2,ZQsL)
 
      S = (LOG(ZTs*USR/VISA/RT)-PTZs)/(von*fdg)
      D = (LOG(ZQs*USR/VISA/RQ)-PQZs)/(von*fdg)
      T_hs =TSR*S +TS -.0098*ZTs
      Q_hs =(QSR*D + QS)*1000
      U_wg_hs = USR*(LOG(ZUs/ZO) - PUZs)/0.4
      IF(U_wg_hs.GE.Wg) THEN
         U_hs = SQRT(U_wg_hs**2 - Wg**2)
      ELSE
         U_hs = U_wg_hs
      ENDIF
c
c Alternatively, you could add the delta correction to the top measurement.
c It shouldn't make a difference as long as log profiles are forced to go
c through both measurements.
c       CALL ZETA(T,Q,USR,TSR,QSR,ZU,ZUL)
c       CALL ZETA(T,Q,USR,TSR,QSR,ZT,ZTL)
c       CALL ZETA(T,Q,USR,TSR,QSR,ZQ,ZQL)
c       PUZ = PSI(1,ZUL)
c       PTZ = PSI(2,ZTL)
c       PQZ = PSI(2,ZQL)
c       U_wg_hs=DU_Wg + USR*(LOG(ZUs/ZU)-(PUZs-PUZ))/0.4
c       U_hs = sqrt(U_wg_hs**2 - Wg**2)
c       T_hs=T-.0098*(ZTs-ZT)+TSR*2.2*(LOG(ZTs/ZT)-PTZs+PTZ)*1.15
c       Q_hs=(Q + QSR*2.2*(LOG(ZQs/ZQ) -PQZs +PQZ)*1.15)*1000
c
      IF(IHUMID.EQ.1) THEN    ! then need to convert sp hum into rh
              Q_hs = Q_hs/1000     ! sh kg/kg
              RHO=1./(287.*(T+273.16)*(1.+.61*Q))*P*100.
              P_hs = P - (RHO*grav*(ZTs - ZT))/100 !Approx hydrost.Pressure mb
              RHO_hs=1./(287.*(T_hs+273.16)*(1.+.61*Q_hs))*P_hs*100
              RHO_avg = (RHO + RHO_hs)/2
              P_hs = P -(RHO_avg*grav*(ZTs - ZT))/100 !hydrostatic Pressure
              call humidity(T_hs,P_hs,QA)         !Teten's formula for Pvap,sat
              ee=Q_hs*P_hs/(.62197 + .378*Q_hs)   !to get vapor pressure
              Q_hs = ee/QA                        !to get relative humidity
      ENDIF
 
      RETURN
      END
c
cc
ccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zang (cosz,rlat,rmonth,riday,hrl)
      integer imd(12)
      data imd/31,28,31,30,31,30,31,31,30,31,30,31/
      save
      month=rmonth
      iday=riday
c      noczon=0
c      fixmon=12.
c       ixday=25
      if (month .gt. 12) month=month-12 
      if (hrl .ge. 24.) then
       hrl=hrl-24.0
       iday=iday+1
         if (iday .gt. imd(month)) then
          iday=iday-imd(month)
          month=month+1
         end if
      end if
       monthi=month
        idayi=iday
c      if (noczon .eq. 1) then
c      monthi=fixmon     
c       idayi= ixday 
c      end if
       fac=.01745
       x=float((monthi-1)*30+idayi+monthi/2)
       xlat=rlat*fac
       phi=2.*3.14159*x/365.
       tphi=phi*2.
       ttphi=phi*3.
       cosphi=cos(phi)
       sinphi=sin(phi)
       costph=cos(tphi)
       sintph=sin(tphi)
       costtp=cos(ttphi)
       sinttp=sin(ttphi)
       delta=.006918-0.399912*cosphi+0.070257*sinphi-0.006758*
     *        costph+0.000907*sintph-0.002697*costtp+.00148*sinttp
       anx=(hrl-12.)*15.*fac
       eq=.000075+.001868*cosphi-.032077*sinphi-.014615*costph
     *     -.040849*sintph
       ang=anx+eq
       cosz=sin(xlat)*sin(delta)+cos(xlat)*cos(delta)*cos(ang)
       rmonth=month
       riday=iday
      write (6,200) monthi,idayi,hrl,cosz
  200 format(1x,2i4,2f10.4)
      return
      end

ccc
cc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RADRAT (IFLAG)
c      IMPLICIT NONE
C     NB IS A NUMBER OF BLOCKS FOR CALCULATION
      INTEGER NBB,NQ,NX,NZ,LAY,NADD,NXI
      INTEGER NW,NX1,NZ2,NZ15,NB1,NB2
      INTEGER IFLAG,NPP1,NY,IMAX,ILES,IL2,JMAX,JLES
      INTEGER JL2,KMAX,KLES,KL2
      PARAMETER (NX=130,NY=130,NZ=43,NBB=1,LAY=88,NADD=7)
      PARAMETER (NQ=NZ+NADD-1,NW=NQ-1)
      PARAMETER (NXI=NX-2,NX1=NXI/NBB)
      PARAMETER (NZ2=NZ*2,NZ15=NZ*15,nb1=nz15+nx,nb2=nz*7)

      REAL    PL(LAY),PA(LAY)

      common/bxyz/ imax,iles,il2,jmax,jles,jl2,kmax,kles,kl2,iw(5)
      common/size/ tnw,tns,tng,roqs,roqg,roqr

      common/sltm/ rlat,rmonth,riday,hrl,so0,cosz,icosz,terman,rsfc
      common/srflx/ sfir(lay),sfne(lay),qswm(lay),qirm(lay)

      common/b5/ tb(nz),qb(nz),zz0(nz),rho(nz),tb1(nz),qb1(nz),zz2(nb1)
      common/b6/ zz3(nz2),p0(nz),pi(nz),zz4(nb2)    
      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      COMMON/SURFACE/ SUN_4(NX,ny,4)

      common/radtemp2/ rsw1(nx,nz),rlw1(nx,nz)
CC
      COMMON/SLWAVE/ RSW(NX,NY,NZ),RLW(NX,NY,NZ)
      REAL  RSIRBM(NX1),RSUVBM(NX1)
      COMMON/CLOUDP/ ICT,ICB
  
C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------
      REAL   TA(LAY),WA(LAY),OA(LAY)
      REAL   PLU(NADD)
      REAL   FCLD(NX1,NW),TAUIR1(NX1,NW),FLX(NX1,NQ)
      REAL   TSQ(NX1,1)
      REAL   PLQ(NX1,NQ),TAQ(NX1,NW),WAQ(NX1,NW)
      REAL   OAQ(NX1,NW)
      REAL   COOLR1(NX1+1,NQ),HEATR1(NX1+1,NQ)
      REAL   REFF(NX1,NW,2)
      REAL   FLX1(NX1,NQ)
      REAL   FDIRIR(NX1),FDIFIR(NX1)
      REAL   FDIRPAR(NX1),FDIFPAR(NX1)
c      LOGICAL HIGH
c      DATA HIGH/.FALSE./
c     data plu/ 0.01                                       /
c     data plu/ 0.01,             4.04,               22.46/
c     data plu/ 0.01      , 3.04      , 7.18        , 22.46/
      data plu/ 0.01, 0.56, 3.04, 4.04, 7.18,  12.38, 22.46/
      INTEGER I,K,IRST,NP,KM,I3,K1,IHALFP,JJ
      REAL    ACO,CMC,SC0
      DATA IRST/0/
      SAVE

      do 2000 jj=2,jles
c
        jm=jj-1
c
      DO K=1,NQ
      DO I=1,nx1
        FLX1(I,K)=0.0
      enddo
      enddo

      DO I=1,NX1
        FDIRIR(I)=0.0
        FDIFIR(I)=0.0
        FDIRPAR(I)=0.0
        FDIFPAR(I)=0.0
      ENDDO


      IF (IRST .NE. 0) GO TO 500
C
C     NP=NUMBER OF ATMOSPHERIC LAYERS; NPP1=SURFACE LEVEL

      NP  =KL2+NADD
      NPP1=NP+1

C  SOLAR CONSTANT AND COSINE OF SOLAR ZENITH ANGLE
C      SC0=1365.
      SC0=1411.

        ICT=24
        ICB=32

        DO I=1,NX1
          RSIRBM(I)=0.07
          RSUVBM(I)=0.07
        ENDDO
C
C  ASSIGN CO2 (CMC). UNITS ARE PARTS/PART
      CMC=300.E-6
      WRITE(6,*) 'NPP1(3)=',NPP1

      DO K=1,KL2+1
        KM=KL2+3-K
        PL(K+NADD)=1.E-3*P00(KM)
      END DO

      DO K=1,NADD
        PL(K)=PLU(K)
      END DO

      DO K=1,KL2
        KM=KL2+2-K
        PA(K+NADD)=1.E-3*P0(KM)
      END DO

      DO 13 K=1,NADD
        PA(K)=0.5*(PL(K)+PL(K+1))
        WA(K)=1.E-6
c        TA(K)=190.
   13 CONTINUE
C     PRINT *,'       RADRAT.F          5  FITO3.F'
      CALL FITO3 (NPP1,PA,TA,OA)
       PRINT*
      PRINT*,' K       PA        PL       AO        TA         WA'
c      WRITE (6,975)(K,PA(K),PL(K),OA(K),TA(K),WA(K),
c     1   K=1,NPP1)
C
      IRST=1
c  975 FORMAT(I3,3F10.3,E12.4,F12.5,F12.7)
  500 CONTINUE

      CALL OPT4 (JM,PL,TA,WA,OA,TAUIR1,fcld,TAQ,WAQ,OAQ,PLQ,TSQ,REFF)

        IF (COSZ .GE. 0.005) THEN
          CALL SORAD (NX1,PLQ,TAQ,WAQ,OAQ,CMC,
     $                REFF,RSIRBM,RSUVBM,COSZ,
     $                FLX1,FDIRIR,FDIFIR,FDIRPAR,FDIFPAR)
        END IF

        CALL IRRAD (NX1,TAUIR1,FCLD,PLQ,TAQ,WAQ,OAQ,CMC,TSQ,FLX)


C SINCE FLX1,ARE NORMALIZED, THEY SHOULD BE MULTIPLIED
C BY SC0*COSZ
C SINCE UPWARD FLUX SHOULD BE POSITIVE, HEATR IS MULTIPLIED BY
C A MINUS SIGN
C SAVE SOLAR AND LONG-WAVE RADIATIVE FLUXES
C
C
        do i=1,il2
           i3=i+1
          sun_4(i3,jj,1)=flx1(i,nw+1)*sc0*cosz
          sun_4(i3,jj,2)=flx(i,nw+1)
        enddo
c
        do k=1,nw
        do i=1,il2
           i3=i+1
          heatr1(i3,k)=(flx1(i,k+1)-flx1(i,k))*8.441874/
     1                   (plq(i,k+1)-plq(i,k))
          heatr1(i3,k)=-heatr1(i3,k)*sc0*cosz
          coolr1(i3,k)=(flx(i,k+1)-flx(i,k))*8.441874/
     1                   (plq(i,k+1)-plq(i,k))
        enddo
        enddo
C
        DO K1=NADD+1,NP
          K=NP+2-K1
        DO I=2,iles
           ACO=1./(86400.*PI(K))
          RSW1(I,K)=ACO*HEATR1(I,K1)
          RLW1(I,K)=ACO*COOLR1(I,K1)
        ENDDO
        ENDDO

         jhalfp=jl2/2+1

      IF (IFLAG .EQ. 2 .and. jj .eq. jhalfp) THEN
        IHALFP=IL2/2+1
        WRITE(6,*) 'COSZ=',COSZ
        WRITE(6,*) 'CHECK RSW1(I,K) AND RLW1(I,K) AT CENTRAL POINT'
        DO K1=NADD+1,NP
          K=NP+2-K1
          WRITE(6,7821) k,pl(k1),RSW1(IHALFP,K),RLW1(IHALFP,K),
     1        heatr1(ihalfp,k1),coolr1(ihalfp,k1)
        ENDDO
      ENDIF

          do 250 i=2,iles
            do k=2,kles
              rsw(i,jj,k)=rsw1(i,k)
              rlw(i,jj,k)=rlw1(i,k)
            enddo
  250     continue

2000  continue


7821  FORMAT(2x,i6,2x,f10.3,2x,4e20.10)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE OPT4(J,PL,TA,WA,OA,TAUIR,FCLD,TAQ,WAQ,OAQ,PLQ,TSQ,REFF)
c      IMPLICIT NONE
C DEFINE VARIABLES AND CALCULATE THE OPTICAL THICKNESS
      INTEGER NX,NZ,NADD,NZ2,NZ15
      INTEGER NB1,NB2,NPP1,NW
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NADD=7)
      PARAMETER (NZ2=NZ*2,NZ15=NZ*15)
      PARAMETER (nb1=nz15+nx,nb2=nz*7,NPP1=NZ+NADD-1)
      PARAMETER (NW=npp1-1)
      PARAMETER (NBB=1)
      PARAMETER (NXI=NX-2,NX1=NXI/NBB)

      REAL  PL(1),TA(1),WA(1),OA(1)
      REAL  TAUIR(NX1,1),REFF(NX1,NW,1)
      COMMON /OPTICAL/ TAUSW(NX1,NW,2)
      REAL  FCLD(NX1,1),TSQ(NX1,1)
      REAL  PLQ(NX1,1),TAQ(NX1,1),WAQ(NX1,1)
      REAL  OAQ(NX1,1)

      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      common/bxyz/ imax,iles,il2,jmax,jles,jl2,kmax,kles,kl2,iw(5)

      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      common/b5/ tb(nz),qb(nz),zz0(nz),rho(nz),tb1(nz),qb1(nz),zz2(nb1)
      common/b6/ zz3(nz2),p0(nz),pi(nz),zz4(nb2)
      COMMON/PICNST/ CPI
C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------
      REAL    AA1(NX,NPP1),BB(NX,NPP1)
      REAL    AA2(NX,NPP1)
      INTEGER I,II,K,KM
      REAL    B2,CPI,EFFRAD
      REAL    Q1,Q2,Q3,Q5
      REAL    TAUQC,TAUQG,TAUQII,TAUQIS,TAUQR,TWCO
      SAVE
      twcz=1.e-6
      tcont=190.
      TWCO=1.E-6
      jp=j+1
      EFFC=0.0015
C
      DO K=1,NPP1
         plk=pl(k)
         DO I=1,il2
            PLQ(I,K)=PLK
         ENDDO
      ENDDO

      DO K=1,NW
         oak=oa(k)
         DO I=1,il2
            OAQ(I,K)=OAK
            FCLD(I,K)=1.0
            REFF(I,K,1)=0.
            REFF(I,K,2)=0.
            TAUIR(I,K)=0.
         ENDDO
      ENDDO

      DO K=1,NADD
         wak=wa(k)
         tak=ta(k)
         DO I=1,il2
            WAQ(I,K)=WAK
            TAQ(I,K)=TAK
         ENDDO
      ENDDO

      DO I=1,il2
         TSQ(I,1)=TAIRSFC(I+1,JP)
      ENDDO

      DO K=2,KLES
         KM=KMAX-K+NADD
         RHODZ0=RHO(K)*DZ0(K)
         qb1k=qb1(k)
         pik=pi(k)
         DO I=1,il2
            II=I+1
            TAUQC=0.0
            TAUQR=0.0
            TAUQIS=0.0
            TAUQII=0.0
            TAUQS=0.0
            TAUQG=0.0
            Q1=QCL(II,JP,K)
            Q2=QRN(II,JP,K)
            Q3=QCI(II,JP,K)
            Q4=QCS(II,JP,K)
            Q5=QCG(II,JP,K)
            WAQ(I,KM)=max(twcz,qb1k+dqv(ii,jp,k))
            TAQ(I,KM)=max(tcont,(tb1(k)+dpt(ii,jp,k))*pik)
            IF(Q1 .GE. TWCO) THEN
               TAUQC=RHODZ0*Q1/EFFC
               REFF(I,KM,2)=EFFC*1.0E4
            ENDIF
            IF(Q2 .GE. TWCO) THEN
               EFFRAD=3./((CPI*TNW*ROQR/(RHO(K)*Q2))**.25)
               TAUQR=RHODZ0*Q2/EFFRAD
            ENDIF
            IF(Q4 .GE. TWCO) THEN
               B3=1.E4*RHO(K)*DZ0(K)*Q4
               EFFRAD=3./((CPI*TNS*ROQS/(RHO(K)*Q4))**.25)
               TAUQS=((RHO(K)*Q4)*DZ0(K))/EFFRAD
            ENDIF
            IF(Q5 .GE. TWCO) THEN
               EFFG=3./((CPI*TNG*ROQG/(RHO(K)*Q5))**.25)
               TAUQG=RHODZ0*Q5/EFFG
            ENDIF
            IF(Q3 .GE. TWCO) THEN
               B2=1.E4*RHO(K)*DZ0(K)*Q3
               EFFRAD=0.0125+(TAQ(I,KM)-243.16)*0.00050
               IF (TAQ(I,KM) .GT. 243.16) EFFRAD=0.0125
               IF (TAQ(I,KM) .LT. 223.16) EFFRAD=0.0025
               TAUQIS=B2*(-0.006656+ 3.686E-4/EFFRAD)
               TAUQII=B2*(-0.011500+ 4.110E-4/EFFRAD
     +               +17.300E-8/(EFFRAD*EFFRAD))
               REFF(I,KM,1)=EFFRAD*1.0E4
            ENDIF
            TAUSW(I,KM,1)=TAUQIS+TAUQG+TAUQS
            TAUSW(I,KM,2)=1.5*(TAUQC+TAUQR)
            TAUIR(I,KM)=0.5*TAUSW(I,KM,2)+TAUQII+TAUQG+TAUQS
         ENDDO
      ENDDO

      km=kmax-2+nadd
      km1=km-1
      DO I=1,il2
         aa1(i,km)=0.8*tausw(i,km,1)+0.2*tausw(i,km1,1)
         aa2(i,km)=0.8*tausw(i,km,2)+0.2*tausw(i,km1,2)
         bb(i,km)=0.8*tauir(i,km)+0.2*tauir(i,km1)
      ENDDO

      DO K=3,KLES-1
         KM=KMAX-K+NADD
         KMM=KM-1
         KMP=KM+1
         DO I=1,il2
            aa1(i,km)=0.2*tausw(i,kmm,1)+0.60*tausw(i,km,1)
     1                 +0.2*tausw(i,kmp,1)
            aa2(i,km)=0.2*tausw(i,kmm,2)+0.60*tausw(i,km,2)
     1                 +0.2*tausw(i,kmp,2)
            bb(i,km)=0.2*tauir(i,kmm)+0.60*tauir(i,km)
     1                +0.2*tauir(i,kmp)
         ENDDO
      ENDDO

      km=kmax-kles+nadd
      km1=km-1
      DO I=1,il2
         aa1(i,km)=0.8*tausw(i,km,1)+0.2*tausw(i,km1,1)
         aa2(i,km)=0.8*tausw(i,km,2)+0.2*tausw(i,km1,2)
         bb(i,km)=0.8*tauir(i,km)+0.2*tauir(i,km1)
      ENDDO

      DO K=2,KLES
         KM=KMAX-K+NADD
         DO I=1,il2
            TAUSW(I,KM,1)=AA1(I,KM)
            TAUSW(I,KM,2)=AA2(I,KM)
            TAUIR(I,KM)=BB(I,KM)
         ENDDO
      ENDDO

      RETURN
      END

C********************* CLIRAD IR1  DATE: OCTOBER, 1994 *************
      SUBROUTINE IRRAD (M,TAUCL,CCLD,PL,TA,WA,OA,CO2,TS,FLX)

      INTEGER NXX,NZ,NADD,NBB,NXI,MM,NP,NX,NO,NT,NPP
      PARAMETER (NX=130,NZ=43,NADD=7,NBB=1)
      PARAMETER (NXI=NX-2)
      PARAMETER (MM=NXI/NBB,NP=NZ-2+NADD,NPP=NP+1)
      PARAMETER (NXX=26,NO=21,NT=7)
      INTEGER M
      REAL    CO2
      REAL    TAUCL(M,1),CCLD(M,1),PL(M,1),
     $        TA(M,1),WA(M,1),OA(M,1),TS(1)
      REAL    FLX(M,1)
      INTEGER IRADAVE,KP,K1P,K1M,K2M
      COMMON/IPTIONR/ IRADAVE
C
C---- STATIC DATA -----

      REAL    CB(5,8)

C---- TEMPORARY ARRAYS -----

      REAL PA(MM,NP),DT(MM,NP)
      REAL SCO3(MM,NPP),SCOPRE(MM,NPP),SCOTEM(MM,NPP)
      REAL DH2O(MM,NP),DCONT(MM,NP),DCO2(MM,NP),DO3(MM,NP)
      REAL TH2O(MM,6),TCON(MM,3),TCO2(MM,6,2)
      REAL H2OEXP(MM,NP,6),CONEXP(MM,NP,3)
      COMMON/radh2o/ H2OEXP
      COMMON/radcon/ CONEXP
      REAL CLR(MM,0:NPP),FCLR(MM)
      REAL BLAYER(MM,0:NPP)
      REAL TRANT(MM)
      REAL FLXU(MM,NPP),FLXD(MM,NPP)

      LOGICAL OZNBND
      LOGICAL CO2BND

      REAL O1 (NXX,NO,NT),O2 (NXX,NO,NT),O3 (NXX,NO,NT)
      REAL DP,XX

       DATA CB/
     1 -2.6844E-1,-8.8994E-2, 1.5676E-3,-2.9349E-6, 2.2233E-9,
     2  3.7315E+1,-7.4758E-1, 4.6151E-3,-6.3260E-6, 3.5647E-9,
     3  3.7187E+1,-3.9085E-1,-6.1072E-4, 1.4534E-5,-1.6863E-8,
     4 -4.1928E+1, 1.0027E+0,-8.5789E-3, 2.9199E-5,-2.5654E-8,
     5 -4.9163E+1, 9.8457E-1,-7.0968E-3, 2.0478E-5,-1.5514E-8,
     6 -1.0345E+2, 1.8636E+0,-1.1753E-2, 2.7864E-5,-1.1998E-8,
     7 -6.9233E+0,-1.5878E-1, 3.9160E-3,-2.4496E-5, 4.9301E-8,
     8  1.1483E+2,-2.2376E+0, 1.6394E-2,-5.3672E-5, 6.6456E-8/

      LOGICAL FIRST
      DATA FIRST /.TRUE./
      INTEGER I,K,IP,IW,IT,IB,IK,IQ,ISB,K1,K2
      INTEGER MMM

      include "o3.tran3"	!cray

      SAVE
      MMM=M

      IF (FIRST) THEN
        DO IW=1,NO      
         DO IP=1,NXX
          O1(IP,IW,1)=1.0-O1(IP,IW,1)
         ENDDO
        ENDDO
        DO IT=2,NT
         DO IW=1,NO
          DO IP=1,NXX
           O1 (IP,IW,IT)= O1(IP,IW,1)
           O2 (IP,IW,IT)= O2(IP,IW,1)
           O3 (IP,IW,IT)= O3(IP,IW,1)
          ENDDO
         ENDDO
        ENDDO
       FIRST=.FALSE.
      ENDIF

      DPC=789.*CO2
      DO K=1,NP
         KP=K+1
         DO I=1,MMM
            PA(I,K)=0.5*(PL(I,K)+PL(I,KP))
            DT(I,K)=TA(I,K)-250.0
            DP           = PL(I,KP)-PL(I,K)
            DH2O (I,K) = 1.02*WA(I,K)*DP+1.E-10
            DCO2 (I,K) = DPC*DP+1.E-10
            DO3  (I,K) = 476.0*OA(I,K)*DP+1.E-10
            XX=PA(I,K)*0.001618*WA(I,K)*WA(I,K)*DP
            DCONT(I,K) = XX*EXP(1800./TA(I,K)-6.081)+1.E-10
            CLR(I,K)=1.0-(CCLD(I,K)*(1.-EXP(-1.66*TAUCL(I,K))))
         ENDDO
      ENDDO

      DO K=1,NPP
        DO I=1,MMM
         FLXU(I,K) = 0.0
         FLXD(I,K) = 0.0
        ENDDO
       ENDDO

c------------------------------------------------------------------------
       DO IB=1,2
       DO K=1,NP
         DO I=1,MMM
          BLAYER(I,K)=TA(I,K)*(TA(I,K)*(TA(I,K)
     *                 *(TA(I,K)*CB(5,IB)+CB(4,IB))+CB(3,IB))
     *                 +CB(2,IB))+CB(1,IB)
         ENDDO
       ENDDO

        DO I=1,MMM
         BLAYER(I,0)   = 0.0
         BLAYER(I,NPP)=TS(I)*(TS(I)*(TS(I)
     *                *(TS(I)*CB(5,IB)+CB(4,IB))+CB(3,IB))
     *                +CB(2,IB))+CB(1,IB)
        ENDDO

      CALL H2OEXPS(IB,DH2O,PA,DT)

      DO 2001 K1=1,NP
         K1P=K1+1
         K1M=K1-1

         DO I=1,MMM
            FCLR(I)=1.0
         ENDDO

        DO IK=1,6
          DO I=1,MMM
           TH2O(I,IK)=1.0
          ENDDO
        ENDDO

      CALL WVKDIS(IB,K1,TCON,TH2O,TRANT)
      DO I=1,MMM
          FCLR(I) = FCLR(I)*CLR(I,K1)
          FLXU(I,K1)=FLXU(I,K1)-BLAYER(I,K1)+(TRANT(I)*(BLAYER(I,K1)
     1              -BLAYER(I,K1P)))*fclr(i)
          FLXD(I,K1P)=FLXD(I,K1P)+BLAYER(I,K1)+(TRANT(I)*(BLAYER(I,K1M)
     1              -BLAYER(I,K1)))*fclr(i)
      ENDDO

      DO 3001 K2=K1P+1,NPP
         K2M=K2-1
         CALL WVKDIS(IB,K2M,TCON,TH2O,TRANT)
         DO I=1,MMM
            FCLR(I) = FCLR(I)*CLR(I,K2M)
            FLXU(I,K1)=FLXU(I,K1)+(TRANT(I)*(BLAYER(I,K2M)
     1                  -BLAYER(I,K2)))*fclr(i)
            FLXD(I,K2)=FLXD(I,K2)+(TRANT(I)*(BLAYER(I,K1M)
     1                  -BLAYER(I,K1)))*fclr(i)
         ENDDO

 3001   CONTINUE
 2001   CONTINUE

       DO I=1,MMM
        FLXU(I,NPP)=FLXU(I,NPP)-BLAYER(I,NPP)
       ENDDO
       ENDDO

c------------------------------------------------------------------------


      DO 1000 IB=3,6
       CO2BND=IB.EQ.3
       OZNBND=IB.EQ.5

       DO K=1,NP
         DO I=1,MMM
          BLAYER(I,K)=TA(I,K)*(TA(I,K)*(TA(I,K)
     *                 *(TA(I,K)*CB(5,IB)+CB(4,IB))+CB(3,IB))
     *                 +CB(2,IB))+CB(1,IB)
         ENDDO
       ENDDO

        DO I=1,MMM
         BLAYER(I,0)   = 0.0
         BLAYER(I,NPP)=TS(I)*(TS(I)*(TS(I)
     *                *(TS(I)*CB(5,IB)+CB(4,IB))+CB(3,IB))
     *                +CB(2,IB))+CB(1,IB)
        ENDDO

      IF (OZNBND) CALL COLUMN(MMM,NP,PA,DT,DO3,SCO3,SCOPRE,SCOTEM)
      CALL H2OEXPS(IB,DH2O,PA,DT)
      CALL CONEXPS(IB,DCONT)
      IF( CO2BND) CALL CO2EXPS(DCO2,PA,DT)

      DO 2000 K1=1,NP
         K1P=K1+1
         K1M=K1-1

         DO I=1,MMM
            FCLR(I)=1.0
         ENDDO

        DO IK=1,6
          DO I=1,MMM
           TH2O(I,IK)=1.0
          ENDDO
        ENDDO

         DO IQ=1,3
           DO I=1,MMM
            TCON(I,IQ)=1.0
           ENDDO
         ENDDO

       IF (CO2BND) THEN
         DO ISB=1,2
          DO IK=1,6
            DO I=1,MMM
             TCO2(I,IK,ISB)=1.0
            ENDDO
          ENDDO
         ENDDO
       ENDIF

      K2M=k1P-1
      CALL WVKDIS(IB,K2M,TCON,TH2O,TRANT)
      IF (CO2BND) CALL CO2KDIS(K2M,TCO2,TRANT)
      IF (OZNBND) CALL TABLUP(K1,K1P,Mmm,NP,NXX,NO,NT,SCO3,SCOPRE,
     *                  SCOTEM,-6.0,-2.0,0.3,0.2,O1,O2,O3,TRANT)
      DO I=1,MMM
          FCLR(I) = FCLR(I)*CLR(I,K2M)
          FLXU(I,K1)=FLXU(I,K1)-BLAYER(I,K1)+(TRANT(I)*(BLAYER(I,K2M)
     1              -BLAYER(I,K1P)))*fclr(i)
          FLXD(I,K1P)=FLXD(I,K1P)+BLAYER(I,K1)+(TRANT(I)*(BLAYER(I,K1M)
     1              -BLAYER(I,K1)))*fclr(i)
      ENDDO

      DO 3000 K2=K1P+1,NPP
         K2M=K2-1
         CALL WVKDIS(IB,K2M,TCON,TH2O,TRANT)
         IF (CO2BND) CALL CO2KDIS(K2M,TCO2,TRANT)
         IF (OZNBND) CALL TABLUP(K1,K2,Mmm,NP,NXX,NO,NT,SCO3,SCOPRE,
     *                  SCOTEM,-6.0,-2.0,0.3,0.2,O1,O2,O3,TRANT)
         DO I=1,MMM
            FCLR(I) = FCLR(I)*CLR(I,K2M)
            FLXU(I,K1)=FLXU(I,K1)+(TRANT(I)*(BLAYER(I,K2M)
     1                  -BLAYER(I,K2)))*fclr(i)
            FLXD(I,K2)=FLXD(I,K2)+(TRANT(I)*(BLAYER(I,K1M)
     1                  -BLAYER(I,K1)))*fclr(i)
         ENDDO

 3000 CONTINUE
 2000 CONTINUE

       DO I=1,MMM
        FLXU(I,NPP)=FLXU(I,NPP)-BLAYER(I,NPP)
       ENDDO

 1000 CONTINUE

c------------------------------------------------------------------------

       DO IB=7,8
       DO K=1,NP
         DO I=1,MMM
          BLAYER(I,K)=TA(I,K)*(TA(I,K)*(TA(I,K)
     *                 *(TA(I,K)*CB(5,IB)+CB(4,IB))+CB(3,IB))
     *                 +CB(2,IB))+CB(1,IB)
         ENDDO
       ENDDO

        DO I=1,MMM
         BLAYER(I,0)   = 0.0
         BLAYER(I,NPP)=TS(I)*(TS(I)*(TS(I)
     *                *(TS(I)*CB(5,IB)+CB(4,IB))+CB(3,IB))
     *                +CB(2,IB))+CB(1,IB)
        ENDDO

      CALL H2OEXPS(IB,DH2O,PA,DT)

      DO 2002 K1=1,NP
         K1P=K1+1
         K1M=K1-1

         DO I=1,MMM
            FCLR(I)=1.0
         ENDDO

        DO IK=1,6
          DO I=1,MMM
           TH2O(I,IK)=1.0
          ENDDO
        ENDDO

      CALL WVKDIS(IB,K1,TCON,TH2O,TRANT)
      DO I=1,MMM
          FCLR(I) = FCLR(I)*CLR(I,K1)
          FLXU(I,K1)=FLXU(I,K1)-BLAYER(I,K1)+(TRANT(I)*(BLAYER(I,K1)
     1              -BLAYER(I,K1P)))*fclr(i)
          FLXD(I,K1P)=FLXD(I,K1P)+BLAYER(I,K1)+(TRANT(I)*(BLAYER(I,K1M)
     1              -BLAYER(I,K1)))*fclr(i)
      ENDDO

      DO 3002 K2=K1P+1,NPP
         K2M=K2-1
         CALL WVKDIS(IB,K2M,TCON,TH2O,TRANT)
         DO I=1,MMM
            FCLR(I) = FCLR(I)*CLR(I,K2M)
            FLXU(I,K1)=FLXU(I,K1)+(TRANT(I)*(BLAYER(I,K2M)
     1                  -BLAYER(I,K2)))*fclr(i)
            FLXD(I,K2)=FLXD(I,K2)+(TRANT(I)*(BLAYER(I,K1M)
     1                  -BLAYER(I,K1)))*fclr(i)
         ENDDO

 3002  CONTINUE
 2002  CONTINUE

       DO I=1,MMM
        FLXU(I,NPP)=FLXU(I,NPP)-BLAYER(I,NPP)
       ENDDO
       ENDDO


      DO K=1,NPP
        DO I=1,MMM
         FLX(I,K)=FLXD(I,K)+FLXU(I,K)
        ENDDO
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C***************  CHARNEY  /SILO/Z1MHY/PENG/SOLAR.F ***10/20/95******
      SUBROUTINE SORAD (M,PL,TA,WA,OA,CO2,REFF,
     *                 RSIRBM,RSUVBM,COSZ,
     *                 FLX,FDIRIR,FDIFIR,FDIRPAR,FDIFPAR)

c      IMPLICIT NONE

      INTEGER NXX,NZ,NADD,NP,NXI,NBB,MM
      PARAMETER (NXX=130,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER (NXI=NXX-2)
      PARAMETER (NBB=1,MM=NXI/NBB)

C-----INPUT PARAMETERS

      INTEGER M
      REAL    PL(M,NP+1),TA(M,NP),WA(M,NP),OA(M,NP)
      REAL    REFF(M,NP,2)
      REAL    RSIRBM(1),RSUVBM(1),CO2
      REAL    TAUCLD(MM,NP,2)
      COMMON /OPTICAL/ TAUCLD

C-----OUTPUT PARAMETERS

      REAL    FLX(M,NP+1), cosz
      REAL    FDIRIR(1),FDIFIR(1)
      REAL    FDIRPAR(1),FDIFPAR(1)

 
      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C-----TEMPORARY ARRAY
 
      INTEGER I,K,MMM,kp
      REAL    DP(MM,NP),WH(MM,NP),OH(MM,NP),SCAL(MM,NP)
      REAL    SWH(MM,NP+1),SO2(MM,NP+1),DF(MM,NP+1)
      REAL    SDF(MM),SCLR(MM),CSM(MM),TAUX,X,EXPMN
      real csmc

      SAVE
 
C-----------------------------------------------------------------
      MMM=M

      CSMC=35./SQRT(1224.*COSZ*COSZ+1.)
       DO I= 1, MMM 
         SWH(I,1)=0. 
         SO2(I,1)=0. 
         CSM(I)=CSMC
      ENDDO

      DO K= 1,NP
         KP=K+1
         DO I= 1, MMM
          DP(I,K)=PL(I,KP)-PL(I,K)
          SCAL(I,K)=DP(I,K)*(.5*(PL(I,K)+PL(I,KP))/300.)**.8
          WH(I,K)=1.02*WA(I,K)*SCAL(I,K)*(1.+0.00135*(TA(I,K)-240.))
          OH(I,K)=1.02*OA(I,K)*DP(I,K)*466.7
          SWH(I,KP)=SWH(I,K)+WH(I,K)
          DF(I,K)=0.
        ENDDO
      ENDDO

        DO I=1, MMM
          DF(I,NP+1)=0.
        ENDDO

      CALL SOLIR (MMM,WH,TAUCLD,REFF,CSM,RSIRBM,FLX,FDIRIR,FDIFIR)
      CALL SOLUV (MMM,OH,DP,TAUCLD,REFF,CSM,RSUVBM,FLX,FDIRPAR,FDIFPAR)

      DO K= 1, NP
        KP=K+1
        DO I= 1, MMM
          SO2(I,KP)=SO2(I,K)+165.22*SCAL(I,K)
          X=SO2(I,KP)*CSM(I)
          DF(I,KP)=DF(I,KP)+0.0287*(1.-EXPMN(-0.00027*SQRT(X)))
        ENDDO
      ENDDO

      DO K= 1, NP
       KP=K+1
        DO I= 1, MMM
         SO2(I,KP)=SO2(I,K)+CO2*789.*SCAL(I,K)
        ENDDO
      ENDDO

      CALL FLXCO2(MMM,NP,SO2,SWH,CSM,DF)

       DO I=1,MMM
          SDF(I)=0.0 
          SCLR(I)=1.0 
       ENDDO

      DO K=1,NP
        KP=K+1
        DO I=1,MMM
           TAUX=TAUCLD(I,K,1)+TAUCLD(I,K,2)
           IF(TAUX.GT.0.01 .AND. SCLR(I).EQ.1.) THEN
              SDF(I)=DF(I,K)
              SCLR(I)=0.0
           ENDIF
           FLX(I,KP)=FLX(I,KP)-SDF(I)-DF(I,KP)*SCLR(I)
        ENDDO
      ENDDO

       DO I= 1, MMM
           FDIRIR(I)=FDIRIR(I)-SDF(I)-DF(I,NP+1)*SCLR(I)
       ENDDO

      RETURN
      END  


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COLUMN (M,NP,PA,DT,SABS0,SABS,SPRE,STEM)
C
C**************************************************************************
C-----COMPUTE COLUMN-INTEGRATED (FROM TOP OF THE MODEL ATMOSPHERE)
C     ABSORBER AMOUNT (SABS), ABSORBER-WEIGHTED PRESSURE (SPRE) AND
C     TEMPERATURE (STEM).
C     COMPUTATIONS OF SPRE AND STEM FOLLOWS EQS. (37) AND (38).
C
C--- INPUT PARAMETERS
C   NUMBER OF SOUNDINGS IN ZONAL DIRECTION (M)
C   NUMBER OF SOUNDINGS IN MERIDIONAL DIRECTION (N)
C   NUMBER OF ATMOSPHERIC LAYERS (NP)
C   LAYER PRESSURE (PA)
C   LAYER TEMPERATURE MINUS 250K (DT)
C   LAYER ABSORBER AMOUNT (SABS0)
C
C--- OUTPUT PARAMETERS
C   COLUMN-INTEGRATED ABSORBER AMOUNT (SABS)
C   COLUMN ABSORBER-WEIGHTED PRESSURE (SPRE)
C   COLUMN ABSORBER-WEIGHTED TEMPERATURE (STEM)
C
C--- UNITS OF PA AND DT ARE MB AND K, RESPECTIVELY.
C    UNITS OF SABS ARE G/CM**2 FOR WATER VAPOR AND (CM-ATM)STP FOR CO2 AND O3
C**************************************************************************
c      IMPLICIT NONE
      INTEGER M,NP

C---- INPUT PARAMETERS -----

      REAL    PA(M,NP),DT(M,NP),SABS0(M,NP)

C---- OUTPUT PARAMETERS -----

      REAL    SABS(M,NP+1),SPRE(M,NP+1),STEM(M,NP+1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

      INTEGER I,K,MMM

      SAVE

C*********************************************************************
CTAO
      MMM=M
         DO I=1,MMM
          SABS(I,1)=0.0
          SPRE(I,1)=0.0
          STEM(I,1)=0.0
        ENDDO

        DO K=1,NP
          KP=K+1
          DO I=1,MMM
           SABS(I,KP)=SABS(I,K)+SABS0(I,K)
           SPRE(I,KP)=SPRE(I,K)+PA(I,K)*SABS0(I,K)
           STEM(I,KP)=STEM(I,K)+DT(I,K)*SABS0(I,K)
          ENDDO
        ENDDO

       RETURN
       END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TABLUP(K1,K2,M,NP,NX,NH,NT,SABS,SPRE,STEM,W1,P1,
     *                  DWE,DPE,COEF1,COEF2,COEF3,TRAN)
C**********************************************************************
C   COMPUTE WATER VAPOR, CO2, AND O3 TRANSMITTANCES BETWEEN LEVELS K1 AND K2
C   USING TABLE LOOK-UP FOR M X N SOUNDINGS.
C
C   CALCULATIONS FOLLOW EQ. (40) OF CHOU AND SUAREZ (1995)
C
C---- INPUT ---------------------
C  INDICES FOR PRESSURE LEVELS (K1 AND K2)
C  NUMBER OF GRID INTERVALS IN ZONAL DIRECTION (M)
C  NUMBER OF GRID INTERVALS IN MERIDIONAL DIRECTION (N)
C  NUMBER OF ATMOSPHERIC LAYERS (NP)
C  NUMBER OF PRESSURE INTERVALS IN THE TABLE (NX)
C  NUMBER OF ABSORBER AMOUNT INTERVALS IN THE TABLE (NH)
C  NUMBER OF TABLES COPIED (NT)
C  COLUMN-INTEGRATED ABSORBER AMOUNT (SABS)
C  COLUMN ABSORBER AMOUNT-WEIGHTED PRESSURE (SPRE)
C  COLUMN ABSORBER AMOUNT-WEIGHTED TEMPERATURE (STEM)
C  FIRST VALUE OF ABSORBER AMOUNT (LOG10) IN THE TABLE (W1) 
C  FIRST VALUE OF PRESSURE (LOG10) IN THE TABLE (P1) 
C  SIZE OF THE INTERVAL OF ABSORBER AMOUNT (LOG10) IN THE TABLE (DWE)
C  SIZE OF THE INTERVAL OF PRESSURE (LOG10) IN THE TABLE (DPE)
C  PRE-COMPUTED COEFFICIENTS (COEF1, COEF2, AND COEF3)
C
C---- UPDATED ---------------------
C  TRANSMITTANCE (TRAN)
C
C  NOTE:
C   (1) UNITS OF SABS ARE G/CM**2 FOR WATER VAPOR AND (CM-ATM)STP FOR CO2 AND O3.
C   (2) UNITS OF SPRE AND STEM ARE, RESPECTIVELY, MB AND K.
C   (3) THERE ARE NT IDENTICAL COPIES OF THE TABLES (COEF1, COEF2, AND
C       COEF3).  THE PRUPOSE OF USING THE MULTIPLE COPIES OF TABLES IS
C       TO INCREASE THE SPEED IN PARALLEL (VECTORIZED) COMPUTATIONS.
C       IF SUCH ADVANTAGE DOES NOT EXIST, NT CAN BE SET TO 1.
C   
C**********************************************************************
c      IMPLICIT NONE
      INTEGER K1,K2,M,NP,NX,NH,NT

C---- INPUT PARAMETERS -----

      REAL    W1,P1,DWE,DPE
      REAL    SABS(M,NP+1),SPRE(M,NP+1),STEM(M,NP+1)
      REAL    COEF1(NX,NH,NT),COEF2(NX,NH,NT),COEF3(NX,NH,NT)

C---- UPDATE PARAMETER -----

      REAL TRAN(1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C---- TEMPORARY VARIABLES -----

      INTEGER I,MMM,IW,IP,NN
      REAL    X1,X2,X3,WE,PE,FW,FP,PA,PB,PC,AX,BA,BB,T1,CA,CB,T2

      SAVE

C**********************************************************************
CTAO
      MMM=M
       DO I=1,MMM

        NN=MOD(I,NT)+1

        X1=SABS(I,K2)-SABS(I,K1)
        X2=(SPRE(I,K2)-SPRE(I,K1))/X1
        X3=(STEM(I,K2)-STEM(I,K1))/X1

        WE=(LOG10(X1)-W1)/DWE
        PE=(LOG10(X2)-P1)/DPE

        WE=MAX(WE,W1-2.*DWE)
        PE=MAX(PE,P1)

        IW=WE+1.5
        IP=PE+1.5

        IW=MIN(IW,NH-1)
        IW=MAX(IW, 2)
        IWM=IW-1
        IWP=IW+1

        IP=MIN(IP,NX-1)
        IP=MAX(IP, 1)
        IPP=IP+1

        FW=WE-IWM
        FP=PE-(IP-1)
        FPM=1.-FP

        PA = COEF1(IP,IWM,NN)*(FPM)+COEF1(IPP,IWM,NN)*FP
        PB = COEF1(IP,IW,  NN)*(FPM)+COEF1(IPP,IW,  NN)*FP
        PC = COEF1(IP,IWP,NN)*(FPM)+COEF1(IPP,IWP,NN)*FP
        AX = (-PA*(1.-FW)+PC*(1.+FW)) *FW*0.5 + PB*(1.-FW*FW)
        BA = COEF2(IP,IW,  NN)*(FPM)+COEF2(IPP,IW,  NN)*FP
        BB = COEF2(IP,IWP,NN)*(FPM)+COEF2(IPP,IWP,NN)*FP
        T1 = BA*(1.-FW) + BB*FW
        CA = COEF3(IP,IW,  NN)*(FPM)+COEF3(IPP,IW,  NN)*FP
        CB = COEF3(IP,IWP,NN)*(FPM)+COEF3(IPP,IWP,NN)*FP
        T2 = CA*(1.-FW) + CB*FW
        TRAN(I)=(AX+(T1+T2*X3)*X3)*TRAN(I)

      ENDDO

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE WVKDIS(IB,K,TCON,TH2O,TRAN)
C**********************************************************************
C   COMPUTE WATER VAPOR TRANSMITTANCE BETWEEN LEVELS K1 AND K2 FOR
C   M X N SOUNDINGS USING THE K-DISTRIBUTION METHOD.
C
C   COMPUTATIONS FOLLOW EQS. (34), (46), (50) AND (52).
C
C---- INPUT PARAMETERS
C  SPECTRAL BAND (IB)
C  NUMBER OF GRID INTERVALS IN ZONAL DIRECTION (M)
C  NUMBER OF GRID INTERVALS IN MERIDIONAL DIRECTION (N)
C  NUMBER OF LEVELS (NP)
C  CURRENT LEVEL (K)
C  EXPONENTIALS FOR LINE ABSORPTION (H2OEXP) 
C  EXPONENTIALS FOR CONTINUUM ABSORPTION (CONEXP) 
C
C---- UPDATED PARAMETERS
C  TRANSMITTANCE BETWEEN LEVELS K1 AND K2 DUE TO
C    WATER VAPOR LINE ABSORPTION (TH2O)
C  TRANSMITTANCE BETWEEN LEVELS K1 AND K2 DUE TO
C    WATER VAPOR CONTINUUM ABSORPTION (TCON)
C  TOTAL TRANSMITTANCE (TRAN)
C
C**********************************************************************
c      IMPLICIT NONE
      INTEGER NXX,NZ,NADD,NBB,NXI,IB,M,NP,K

      PARAMETER (NXX=130,NZ=43,NADD=7,NBB=1)
      PARAMETER (NXI=NXX-2)
      PARAMETER (M=NXI/NBB,NP=NZ-2+NADD)

C---- INPUT PARAMETERS ------

      REAL    CONEXP(M,NP,3)
      REAL    H2OEXP(M,NP,6)
      common/radh2o/ H2OEXP
      common/radcon/ CONEXP

C---- UPDATED PARAMETERS -----

      REAL    TH2O(M,1),TCON(M,1),TRAN(1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C---- STATIC DATA -----
C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------

      INTEGER NE(8)
      REAL    FKW(6,8),GKW(6,3)


C-----FKW IS THE PLANCK-WEIGHTED K-DISTRIBUTION FUNCTION DUE TO H2O
C     LINE ABSORPTION GIVEN IN TABLE 4 OF CHOU AND SUAREZ (1995).
C     THE K-DISTRIBUTION FUNCTION FOR THE THIRD BAND, FKW(*,3), IS NOT USED
 
      DATA FKW / 0.2747,0.2717,0.2752,0.1177,0.0352,0.0255,
     2           0.1521,0.3974,0.1778,0.1826,0.0374,0.0527,
     3           6*1.00,
     4           0.4654,0.2991,0.1343,0.0646,0.0226,0.0140,
     5           0.5543,0.2723,0.1131,0.0443,0.0160,0.0000,
     6           0.1846,0.2732,0.2353,0.1613,0.1146,0.0310,
     7           0.0740,0.1636,0.4174,0.1783,0.1101,0.0566,
     8           0.1437,0.2197,0.3185,0.2351,0.0647,0.0183/

C-----GKW IS THE PLANCK-WEIGHTED K-DISTRIBUTION FUNCTION DUE TO H2O
C     LINE ABSORPTION IN THE 3 SUBBANDS (800-720,620-720,540-620 /CM)
C     OF BAND 3 GIVEN IN TABLE 7.  NOTE THAT THE ORDER OF THE SUB-BANDS
C     IS REVERSED.

      DATA GKW/  0.1782,0.0593,0.0215,0.0068,0.0022,0.0000,
     2           0.0923,0.1675,0.0923,0.0187,0.0178,0.0000,
     3           0.0000,0.1083,0.1581,0.0455,0.0274,0.0041/

 
 
C-----NE IS THE NUMBER OF TERMS USED IN EACH BAND TO COMPUTE WATER VAPOR
C     CONTINUUM TRANSMITTANCE (TABLE 6).
 
      DATA NE /0,0,3,1,1,1,0,0/


C-----TCO2 ARE THE SIX EXP FACTORS BETWEEN LEVELS K1 AND K2 
C     TRAN IS THE UPDATED TOTAL TRANSMITTANCE BETWEEN LEVELS K1 AND K2


C-----TH2O IS THE 6 EXP FACTORS BETWEEN LEVELS K1 AND K2 DUE TO
C     H2O LINE ABSORPTION. 

C-----TCON IS THE 3 EXP FACTORS BETWEEN LEVELS K1 AND K2 DUE TO
C     H2O CONTINUUM ABSORPTION.

C-----TRNTH2O IS THE TOTAL TRANSMITTANCE BETWEEN LEVELS K1 AND K2 DUE
C     TO BOTH LINE AND CONTINUUM ABSORPTION COMPUTED FROM EQ. (52).


      INTEGER I,MMM

      SAVE

      MMM=M

      IF (NE(IB).EQ.0) THEN

          DO I=1,MMM

           TH2O(I,1) = TH2O(I,1)*H2OEXP(I,K,1)
           TH2O(I,2) = TH2O(I,2)*H2OEXP(I,K,2)
           TH2O(I,3) = TH2O(I,3)*H2OEXP(I,K,3)
           TH2O(I,4) = TH2O(I,4)*H2OEXP(I,K,4)
           TH2O(I,5) = TH2O(I,5)*H2OEXP(I,K,5)
           TH2O(I,6) = TH2O(I,6)*H2OEXP(I,K,6)

           TRAN(I)      =(FKW(1,IB)*TH2O(I,1)
     *                  + FKW(2,IB)*TH2O(I,2)
     *                  + FKW(3,IB)*TH2O(I,3)
     *                  + FKW(4,IB)*TH2O(I,4)
     *                  + FKW(5,IB)*TH2O(I,5)
     *                  + FKW(6,IB)*TH2O(I,6))
          ENDDO


      ELSEIF (NE(IB).EQ.1) THEN


          DO I=1,MMM
           TH2O(I,1) = TH2O(I,1)*H2OEXP(I,K,1)
           TH2O(I,2) = TH2O(I,2)*H2OEXP(I,K,2)
           TH2O(I,3) = TH2O(I,3)*H2OEXP(I,K,3)
           TH2O(I,4) = TH2O(I,4)*H2OEXP(I,K,4)
           TH2O(I,5) = TH2O(I,5)*H2OEXP(I,K,5)
           TH2O(I,6) = TH2O(I,6)*H2OEXP(I,K,6)

           TCON(I,1)= TCON(I,1)*CONEXP(I,K,1)

           TRAN(I)      =(FKW(1,IB)*TH2O(I,1)
     *                  + FKW(2,IB)*TH2O(I,2)
     *                  + FKW(3,IB)*TH2O(I,3)
     *                  + FKW(4,IB)*TH2O(I,4)
     *                  + FKW(5,IB)*TH2O(I,5)
     *                  + FKW(6,IB)*TH2O(I,6))*TCON(I,1)
          ENDDO

      ELSE

          DO I=1,MMM

           TH2O(I,1) = TH2O(I,1)*H2OEXP(I,K,1)
           TH2O(I,2) = TH2O(I,2)*H2OEXP(I,K,2)
           TH2O(I,3) = TH2O(I,3)*H2OEXP(I,K,3)
           TH2O(I,4) = TH2O(I,4)*H2OEXP(I,K,4)
           TH2O(I,5) = TH2O(I,5)*H2OEXP(I,K,5)
           TH2O(I,6) = TH2O(I,6)*H2OEXP(I,K,6)

           TCON(I,1)= TCON(I,1)*CONEXP(I,K,1)
           TCON(I,2)= TCON(I,2)*CONEXP(I,K,2)
           TCON(I,3)= TCON(I,3)*CONEXP(I,K,3)

           TRAN(I)    = (  GKW(1,1)*TH2O(I,1)
     *                     + GKW(2,1)*TH2O(I,2)
     *                     + GKW(3,1)*TH2O(I,3)
     *                     + GKW(4,1)*TH2O(I,4)
     *                     + GKW(5,1)*TH2O(I,5) ) * TCON(I,1)
c     *                     + GKW(6,1)*TH2O(I,6) ) * TCON(I,1)
     *                  + (  GKW(1,2)*TH2O(I,1)
     *                     + GKW(2,2)*TH2O(I,2)
     *                     + GKW(3,2)*TH2O(I,3)
     *                     + GKW(4,2)*TH2O(I,4)
     *                     + GKW(5,2)*TH2O(I,5) ) * TCON(I,2)
c     *                     + GKW(6,2)*TH2O(I,6) ) * TCON(I,2)
c     *                  + (  GKW(1,3)*TH2O(I,1)
     *                     + (GKW(2,3)*TH2O(I,2)
     *                     + GKW(3,3)*TH2O(I,3)
     *                     + GKW(4,3)*TH2O(I,4)
     *                     + GKW(5,3)*TH2O(I,5)
     *                     + GKW(6,3)*TH2O(I,6) ) * TCON(I,3)
          ENDDO
      ENDIF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CO2KDIS(K,TCO2,TRAN)
C**********************************************************************
C   COMPUTE CO2 TRANSMITTANCES BETWEEN LEVELS K1 AND K2 FOR M X N SOUNDINGS
C   USING THE K-DISTRIBUTION METHOD WITH LINEAR PRESSURE SCALING.
C
C   COMPUTATIONS FOLLOW EQ. (34).
C
C---- INPUT PARAMETERS
C   NUMBER OF GRID INTERVALS IN ZONAL DIRECTION (M)
C   NUMBER OF GRID INTERVALS IN MERIDIONAL DIRECTION (N)
C
C---- UPDATED PARAMETERS
C   TRANSMITTANCE BETWEEN LEVELS K1 AND K2 DUE TO CO2 ABSORPTION
C     FOR THE VARIOUS VALUES OF THE ABSORPTION COEFFICIENT (TCO2)
C   TOTAL TRANSMITTANCE (TRAN)
C
C**********************************************************************
c      IMPLICIT NONE
      integer nx,nz,nadd,nbb,nxi
      INTEGER M,NP,K

C---- INPUT PARAMETERS -----
      PARAMETER (NX=130,NZ=43,NADD=7,NBB=1)
      PARAMETER (NXI=NX-2)
      PARAMETER (M=NXI/NBB,NP=NZ-2+NADD)

      REAL    CO2EXP(M,NP,6,2)
      common /radco2/ co2exp

C---- UPDATED PARAMETERS -----

      REAL    TCO2(M,6,2),TRAN(1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE
C---- STATIC DATA -----

      REAL    GKC(6,2)

C---- TEMPORARY ARRAYS -----

      REAL    XC

C-----GKC IS THE PLANCK-WEIGHTED CO2 K-DISTRIBUTION FUNCTION 
C     IN THE BAND-WING AND BAND-CENTER REGIONS GIVEN IN TABLE 7.
C     FOR COMPUTING EFFICIENCY, SUB-BANDS 3A AND 3C ARE COMBINED.

      DATA GKC/  0.1395,0.1407,0.1549,0.1357,0.0182,0.0220,
     2           0.0766,0.1372,0.1189,0.0335,0.0169,0.0059/

C-----TCO2 IS THE 6 EXP FACTORS BETWEEN LEVELS K1 AND K2. 
C     XC IS THE TOTAL CO2 TRANSMITTANCE GIVEN BY EQ. (53).

      INTEGER I,MMM

      SAVE

CTAO
      MMM=M
          DO I=1,MMM

C-----BAND-WINGS

           TCO2(I,1,1)=TCO2(I,1,1)*CO2EXP(I,K,1,1)
           XC=             GKC(1,1)*TCO2(I,1,1)

           TCO2(I,2,1)=TCO2(I,2,1)*CO2EXP(I,K,2,1)
           XC=XC+GKC(2,1)*TCO2(I,2,1)

           TCO2(I,3,1)=TCO2(I,3,1)*CO2EXP(I,K,3,1)
           XC=XC+GKC(3,1)*TCO2(I,3,1)

           TCO2(I,4,1)=TCO2(I,4,1)*CO2EXP(I,K,4,1)
           XC=XC+GKC(4,1)*TCO2(I,4,1)

           TCO2(I,5,1)=TCO2(I,5,1)*CO2EXP(I,K,5,1)
           XC=XC+GKC(5,1)*TCO2(I,5,1)

           TCO2(I,6,1)=TCO2(I,6,1)*CO2EXP(I,K,6,1)
           XC=XC+GKC(6,1)*TCO2(I,6,1)

C-----BAND-CENTER REGION

           TCO2(I,1,2)=TCO2(I,1,2)*CO2EXP(I,K,1,2)
           XC=XC+GKC(1,2)*TCO2(I,1,2)

           TCO2(I,2,2)=TCO2(I,2,2)*CO2EXP(I,K,2,2)
           XC=XC+GKC(2,2)*TCO2(I,2,2)

           TCO2(I,3,2)=TCO2(I,3,2)*CO2EXP(I,K,3,2)
           XC=XC+GKC(3,2)*TCO2(I,3,2)

           TCO2(I,4,2)=TCO2(I,4,2)*CO2EXP(I,K,4,2)
           XC=XC+GKC(4,2)*TCO2(I,4,2)

           TCO2(I,5,2)=TCO2(I,5,2)*CO2EXP(I,K,5,2)
           XC=XC+GKC(5,2)*TCO2(I,5,2)

           TCO2(I,6,2)=TCO2(I,6,2)*CO2EXP(I,K,6,2)
           XC=XC+GKC(6,2)*TCO2(I,6,2)

           TRAN(I)=TRAN(I)*XC

         ENDDO

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE H2OEXPS(IB,DH2O,PA,DT)
C**********************************************************************
C   COMPUTE EXPONENTIALS FOR WATER VAPOR LINE ABSORPTION
C   IN INDIVIDUAL LAYERS.
C
C---- INPUT PARAMETERS
C  SPECTRAL BAND (IB)
C  NUMBER OF GRID INTERVALS IN ZONAL DIRECTION (M)
C  NUMBER OF GRID INTERVALS IN MERIDIONAL DIRECTION (N)
C  NUMBER OF LAYERS (NP)
C  LAYER WATER VAPOR AMOUNT FOR LINE ABSORPTION (DH2O) 
C  LAYER PRESSURE (PA)
C  LAYER TEMPERATURE MINUS 250K (DT)
C
C---- OUTPUT PARAMETERS
C  6 EXPONENTIALS FOR EACH LAYER  (H2OEXP)
C
C**********************************************************************
c      IMPLICIT NONE
      INTEGER NXX,NZ,NADD,NBB,NXI,IB,M,NP
      PARAMETER (NXX=130,NZ=43,NADD=7,NBB=1)
      PARAMETER (NXI=NXX-2)
      PARAMETER (M=NXI/NBB,NP=NZ-2+NADD)

C---- INPUT PARAMETERS ------

      REAL    DH2O(M,NP),PA(M,NP),DT(M,NP)

C---- OUTPUT PARAMETERS -----

      REAL    H2OEXP(M,NP,6)
      COMMON/radh2o/ H2OEXP

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C---- STATIC DATA -----

      INTEGER MW(8)
      REAL    XKW(8),AW(8),BW(8)

C---- TEMPORARY ARRAYS -----

      REAL    XH

C-----XKW  ARE THE ABSORPTION COEFFICIENTS FOR THE FIRST
C     K-DISTRIBUTION FUNCTION DUE TO WATER VAPOR LINE ABSORPTION
C     (TABLES 4 AND 7).  UNITS ARE CM**2/G    
 
      DATA XKW / 29.55  , 4.167E-1, 1.328E-2, 5.250E-4,
     *            5.25E-4, 2.340E-3, 1.320E-0, 5.250E-4/
 
C-----MW ARE THE RATIOS BETWEEN NEIGHBORING ABSORPTION COEFFICIENTS
C     FOR WATER VAPOR LINE ABSORPTION (TABLES 4 AND 7).
 
      DATA MW /6,6,8,6,6,8,6,16/

C-----AW AND BW (TABLE 3) ARE THE COEFFICIENTS FOR TEMPERATURE SCALING
C     IN EQ. (25).
 
      DATA AW/ 0.0021, 0.0140, 0.0167, 0.0302,
     *         0.0307, 0.0154, 0.0008, 0.0096/
      DATA BW/ -1.01E-5, 5.57E-5, 8.54E-5, 2.96E-4,
     *          2.86E-4, 7.53E-5,-3.52E-6, 1.64E-5/

      INTEGER I,K,IK,MMM

      SAVE

C**********************************************************************
C    NOTE THAT THE 3 SUB-BANDS IN BAND 3 USE THE SAME SET OF XKW, AW,
C    AND BW.  THEREFORE, H2OEXP FOR THESE SUB-BANDS ARE IDENTICAL.
C**********************************************************************
CTAO
      MMM=M
 
        DO K=1,NP
          DO I=1,MMM
           XH = DH2O(I,K)*(PA(I,K)*0.002)
     1        * ( 1.+(AW(IB)+BW(IB)* DT(I,K))*DT(I,K) )
           H2OEXP(I,K,1) = EXP(-XH*XKW(IB))
          ENDDO
        ENDDO

        DO IK=2,6
         IKM=IK-1

         IF(MW(IB).EQ.6) THEN

          DO K=1,NP
            DO I=1,MMM
             XH = H2OEXP(I,K,IKM)*H2OEXP(I,K,IKM)
             H2OEXP(I,K,IK) = XH*XH*XH
            ENDDO
          ENDDO

        ELSEIF(MW(IB).EQ.8) THEN

          DO K=1,NP
            DO I=1,MMM
             XH = H2OEXP(I,K,IKM)*H2OEXP(I,K,IKM)
             XH = XH*XH
             H2OEXP(I,K,IK) = XH*XH
            ENDDO
          ENDDO

        ELSE

          DO K=1,NP
            DO I=1,MMM
             XH = H2OEXP(I,K,IKM)*H2OEXP(I,K,IKM)
             XH = XH*XH
             XH = XH*XH
             H2OEXP(I,K,IK) = XH*XH
            ENDDO
          ENDDO

        ENDIF
       ENDDO


      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CONEXPS(IB,DCONT)
C**********************************************************************
C   COMPUTE EXPONENTIALS FOR CONTINUUM ABSORPTION IN INDIVIDUAL LAYERS.
C
C---- INPUT PARAMETERS
C  SPECTRAL BAND (IB)
C  NUMBER OF GRID INTERVALS IN ZONAL DIRECTION (M)
C  NUMBER OF GRID INTERVALS IN MERIDIONAL DIRECTION (N)
C  NUMBER OF LAYERS (NP)
C  LAYER SCALED WATER VAPOR AMOUNT FOR CONTINUUM ABSORPTION (DCONT) 
C
C---- OUTPUT PARAMETERS
C  1 OR 3 EXPONENTIALS FOR EACH LAYER (CONEXP)
C
C**********************************************************************
c      IMPLICIT NONE
      INTEGER NXX,NZ,NADD,NBB,NXI,IB,M,NP
      PARAMETER (NXX=130,NZ=43,NADD=7,NBB=1)
      PARAMETER (NXI=NXX-2)
      PARAMETER (M=NXI/NBB,NP=NZ-2+NADD)

C---- INPUT PARAMETERS ------

      REAL    DCONT(M,NP)

C---- UPDATED PARAMETERS -----

      REAL    CONEXP(M,NP,3)
      COMMON/radcon/ CONEXP

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C---- STATIC DATA -----

      REAL    XKE(8)


C-----XKE ARE THE ABSORPTION COEFFICIENTS FOR THE FIRST
C     K-DISTRIBUTION FUNCTION DUE TO WATER VAPOR CONTINUUM ABSORPTION
C     (TABLE 6).  UNITS ARE CM**2/G
 
      DATA XKE /  0.00,   0.00,   27.40,   15.8,
     *            9.40,   7.75,     0.0,    0.0/
 

      INTEGER I,K,IQ,MMM

      SAVE

C**********************************************************************
CTAO
      MMM=M

        DO K=1,NP
          DO I=1,MMM
           CONEXP(I,K,1) = EXP(-DCONT(I,K)*XKE(IB))
          ENDDO
        ENDDO

       IF (IB .EQ. 3) THEN

C-----THE ABSORPTION COEFFICIENTS FOR SUB-BANDS 3B (IQ=2) AND 3A (IQ=3)
C     ARE, RESPECTIVELY, DOUBLE AND QUADRUPLE THAT FOR SUB-BAND 3C (IQ=1)
C     (TABLE 6).

        DO IQ=2,3
         iqm=iq-1
         DO K=1,NP
           DO I=1,MMM
            CONEXP(I,K,IQ) = CONEXP(I,K,IQM) *CONEXP(I,K,IQM)
           ENDDO
          ENDDO
        ENDDO

       ENDIF

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CO2EXPS(DCO2,PA,DT)
C       
C**********************************************************************
C   COMPUTE CO2 EXPONENTIALS FOR INDIVIDUAL LAYERS.
C
C---- INPUT PARAMETERS
C  NUMBER OF GRID INTERVALS IN ZONAL DIRECTION (M)
C  NUMBER OF GRID INTERVALS IN MERIDIONAL DIRECTION (N)
C  NUMBER OF LAYERS (NP)
C  LAYER CO2 AMOUNT (DCO2)
C  LAYER PRESSURE (PA)
C  LAYER TEMPERATURE MINUS 250K (DT)
C
C---- OUTPUT PARAMETERS
C  6 EXPONENTIALS FOR EACH LAYER (CO2EXP)
C**********************************************************************
c      IMPLICIT NONE
      INTEGER M,NP
      integer nx,nz,nadd,nbb,nxi
      PARAMETER (NX=130,NZ=43,NADD=7,NBB=1)
      PARAMETER (NXI=NX-2)
      PARAMETER (M=NXI/NBB,NP=NZ-2+NADD)

C---- INPUT PARAMETERS -----

      REAL    DCO2(M,NP),PA(M,NP),DT(M,NP)

C---- OUTPUT PARAMETERS -----

      REAL    CO2EXP(M,NP,6,2)
      common /radco2/ co2exp

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C---- STATIC DATA -----

      REAL    XKC(2),AC(2),BC(2),PM(2),PRC(2)

C---- TEMPORARY ARRAYS -----

      REAL    XC

C-----XKC IS THE ABSORPTION COEFFICIENTS FOR THE
C     FIRST K-DISTRIBUTION FUNCTION DUE TO CO2 (TABLE 7).
C     UNITS ARE 1/(CM-ATM)STP.
 
      DATA XKC/2.656E-5,2.656E-3/
 
C-----PARAMETERS (TABLE 3) FOR COMPUTING THE SCALED CO2 AMOUNT
C     USING (27).

      DATA PRC/  300.0,   30.0/
      DATA PM /    0.5,   0.85/
      DATA AC / 0.0182, 0.0042/
      DATA BC /1.07E-4,2.00E-5/
 
      INTEGER I,K,MMM

      SAVE
C**********************************************************************
CTAO
      MMM=M
        DO K=1,NP
          DO I=1,MMM

C-----COMPUTE THE SCALED CO2 AMOUNT FROM EQ. (27) FOR BAND-WINGS
C     (SUB-BANDS 3A AND 3C).

           XC = DCO2(I,K)*(PA(I,K)/PRC(1))**PM(1)
     1             *(1.+(AC(1)+BC(1)*DT(I,K))*DT(I,K))

C-----SIX EXPONENTIAL BY POWERS OF 8 (TABLE 7).

           CO2EXP(I,K,1,1)=EXP(-XC*XKC(1))

           XC=CO2EXP(I,K,1,1)*CO2EXP(I,K,1,1)
           XC=XC*XC
           CO2EXP(I,K,2,1)=XC*XC

           XC=CO2EXP(I,K,2,1)*CO2EXP(I,K,2,1)
           XC=XC*XC
           CO2EXP(I,K,3,1)=XC*XC

           XC=CO2EXP(I,K,3,1)*CO2EXP(I,K,3,1)
           XC=XC*XC
           CO2EXP(I,K,4,1)=XC*XC

           XC=CO2EXP(I,K,4,1)*CO2EXP(I,K,4,1)
           XC=XC*XC
           CO2EXP(I,K,5,1)=XC*XC

           XC=CO2EXP(I,K,5,1)*CO2EXP(I,K,5,1)
           XC=XC*XC
           CO2EXP(I,K,6,1)=XC*XC

C-----COMPUTE THE SCALED CO2 AMOUNT FROM EQ. (27) FOR BAND-CENTER
C     REGION (SUB-BAND 3B).

           XC = DCO2(I,K)*(PA(I,K)/PRC(2))**PM(2)
     1             *(1.+(AC(2)+BC(2)*DT(I,K))*DT(I,K))

           CO2EXP(I,K,1,2)=EXP(-XC*XKC(2))

           XC=CO2EXP(I,K,1,2)*CO2EXP(I,K,1,2)
           XC=XC*XC
           CO2EXP(I,K,2,2)=XC*XC

           XC=CO2EXP(I,K,2,2)*CO2EXP(I,K,2,2)
           XC=XC*XC
           CO2EXP(I,K,3,2)=XC*XC

           XC=CO2EXP(I,K,3,2)*CO2EXP(I,K,3,2)
           XC=XC*XC
           CO2EXP(I,K,4,2)=XC*XC

           XC=CO2EXP(I,K,4,2)*CO2EXP(I,K,4,2)
           XC=XC*XC
           CO2EXP(I,K,5,2)=XC*XC

           XC=CO2EXP(I,K,5,2)*CO2EXP(I,K,5,2)
           XC=XC*XC
           CO2EXP(I,K,6,2)=XC*XC

          ENDDO
        ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLIR (M,WH,TAUCLD,REFF,
     $                  CSM,RSIRBM,FLX,FDIRIR,FDIFIR)

c      IMPLICIT NONE

      INTEGER NXX,NZ,NADD,NP,NXI,NBB,MM,NK,NBAND
      PARAMETER (NXX=130,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER (NXI=NXX-2)
      PARAMETER (NBB=1,MM=NXI/NBB)
      PARAMETER (NK=10,NBAND=3)
C************************************************************************
C  COMPUTE SOLAR FLUX IN THE INFRARED REGION. THE SPECTRUM IS DIVIDED
C   INTO THREE BANDS:
C
C          BAND   WAVENUMBER(/CM)  WAVELENGTH (MICRON)
C           1       1000-4400         2.27-10.0
C           2       4400-8200         1.22-2.27
C           3       8200-14300        0.70-1.22
C
C----- INPUT PARAMETERS:                            UNITS      SIZE
C
C   NUMBER OF SOUNDINGS IN ZONAL DIRECTION (M)       N/D        1
C   NUMBER OF SOUNDINGS IN MERIDIONAL DIRECTION (N)  N/D        1
C   MAXIMUM NUMBER OF SOUNDINGS IN                   N/D        1
C          MERIDIONAL DIRECTION (NDIM)
C   NUMBER OF ATMOSPHERIC LAYERS (NP)                N/D        1
C   LAYER WATER VAPOR CONTENT (WH)                 GM/CM^2    M*N*NP
C   CLOUD OPTICAL THICKNESS (TAUCLD)                 N/D      M*NDIM*NP*2
C          INDEX 1 FOR ICE PATICLES
C          INDEX 2 FOR LIQUID PARTICLES
C   EFFECTIVE CLOUD-PARTICLE SIZE (REFF)           MICROMETER M*NDIM*NP*2
C          INDEX 1 FOR ICE PATICLES
C          INDEX 2 FOR LIQUID PARTICLES
C   AEROSOL OPTICAL THICKNESS (TAUAL)                N/D      M*NDIM*NP
C   COSECANT OF THE SOLAR ZENITH ANGLE (CSM)         N/D      M*N
C   NEAR IR SURFACE ALBEDO FOR BEAM                FRACTION   M*NDIM
C                RADIATION (RSIRBM)
C
C----- OUTPUT (UPDATED) PARAMETERS:
C
C   ALL-SKY FLUX (DOWNWARD-UPWARD) (FLX)           FRACTION   M*NDIM*(NP+1)
C   ALL-SKY DIRECT DOWNWARD IR FLUX AT
C          THE SURFACE (FDIRIR)                    FRACTION   M*NDIM
C   ALL-SKY DIFFUSE DOWNWARD IR FLUX AT
C          THE SURFACE (FDIFIR)                    FRACTION   M*NDIM
C
C----- NOTE: THE FOLLOWING PARAMETERS MUST BE SPECIFIED BY USERS:
C   AEROSOL SINGLE SCATTERING ALBEDO (SSAAL)         N/D      NBAND
C   AEROSOL ASYMMETRY FACTOR (ASYAL)                 N/D      NBAND
C
C*************************************************************************
C-----INPUT PARAMETERS

      INTEGER M
      REAL    TAUCLD(M,NP,2),REFF(M,NP,2)
      REAL    RSIRBM(1)
      REAL    WH(M,NP),CSM(1)

C-----OUTPUT (UPDATED) PARAMETERS

      REAL    FLX(M,NP+1)
      REAL    FDIRIR(1),FDIFIR(1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C-----STATIC PARAMETERS

      REAL    XK(NK),HK(NBAND,NK)
      REAL    AIA(NBAND,3),AWA(NBAND,3),AIG(NBAND,3),AWG(NBAND,3)

C-----TEMPORARY ARRAY

      INTEGER IB,IK,I,K,MMM
      REAL    SSACL(MM,NP),ASYCL(MM,NP)
      REAL    RR(MM,NP+1),TT(MM,NP+1),TD(MM,NP+1),
     $        RS(MM,NP+1),TS(MM,NP+1)
      REAL    FLXDN(MM,1,NP+1),FDNDIR(MM),FDNDIF(MM)
      REAL    TAUSTO,TAUTO,SSATO,ASYTO
      REAL    TAUX(MM,NP),REFF1,REFF2,W1,W2,G1,G2
 
      real ff,zth,xx,st1,st4,st2,ell,fll,all,bll,cll,dll
      real alf2,alf1,gm3,st7,st3,gm1,gm2,taup,sscp,gp,akk,st8

C-----WATER VAPOR ABSORPTION COEFFICIENT FOR 10 K-INTERVALS.
C     UNIT: CM^2/GM

      DATA XK/            
     1  0.0010, 0.0133, 0.0422, 0.1334, 0.4217,            
     2  1.334,  5.623,  31.62,  177.8,  1000.0/  

C-----WATER VAPOR K-DISTRIBUTION FUNCTION,
C     THE SUM OF HK IS 0.52926. UNIT: FRACTION

      DATA HK/
     1 .01074,.08236,.20673,  .00360,.01157,.03497,
     2 .00411,.01133,.03011,  .00421,.01143,.02260,
     3 .00389,.01240,.01336,  .00326,.01258,.00696,
     4 .00499,.01381,.00441,  .00465,.00650,.00115,
     5 .00245,.00244,.00026,  .00145,.00094,.00000/

 
C-----COEFFICIENTS FOR COMPUTING THE SINGLE SCATTERING ALBEDO OF
C     ICE CLOUDS FROM SSA=1-(AIA(*,1)+AIA(*,2)*REFF+AIA(*,3)*REFF**2)

      DATA AIA/
     1  .08938331, .00215346,-.00000260,
     2  .00299387, .00073709, .00000746,
     3 -.00001038,-.00000134, .00000000/

C-----COEFFICIENTS FOR COMPUTING THE SINGLE SCATTERING ALBEDO OF
C     LIQUID CLOUDS FROM SSA=1-(AWA(*,1)+AWA(*,2)*REFF+AWA(*,3)*REFF**2)

      DATA AWA/
     1  .01209318,-.00019934, .00000007,
     2  .01784739, .00088757, .00000845,
     3 -.00036910,-.00000650,-.00000004/

C-----COEFFICIENTS FOR COMPUTING THE ASYMMETRY FACTOR OF ICE CLOUDS
C     FROM ASYCL=AIG(*,1)+AIG(*,2)*REFF+AIG(*,3)*REFF**2

      DATA AIG/
     1  .84090400, .76098937, .74935228,
     2  .00126222, .00141864, .00119715,
     3 -.00000385,-.00000396,-.00000367/

C-----COEFFICIENTS FOR COMPUTING THE ASYMMETRY FACTOR OF LIQUID CLOUDS
C     FROM ASYCL=AWG(*,1)+AWG(*,2)*REFF+AWG(*,3)*REFF**2

      DATA AWG/
     1  .83530748, .74513197, .79375035,
     2  .00257181, .01370071, .00832441,
     3  .00005519,-.00038203,-.00023263/

      SAVE

C-----INITIALIZE SURFACE FLUXES, REFLECTANCES, AND TRANSMITTANCES
CTAO
      MMM=M
CTAO

C-----INTEGRATION OVER SPECTRAL BANDS

      DO 100 IB=1,NBAND

C-----COMPUTE CLOUD SINGLE SCATTERING ALBEDO AND ASYMMETRY FACTOR
C     FOR A MIXTURE OF ICE AND LIQUID PARTICLES.

       AIA1=AIA(IB,1)
       AIA2=AIA(IB,2)
       AIA3=AIA(IB,3)
       AWA1=AWA(IB,1)
       AWA2=AWA(IB,2)
       AWA3=AWA(IB,3)
       AIG1=AIG(IB,1)
       AIG2=AIG(IB,2)
       AIG3=AIG(IB,3)
       AWG1=AWG(IB,1)
       AWG2=AWG(IB,2)
       AWG3=AWG(IB,3)
       DO K= 1, NP
         DO I= 1, MMM
           SSACL(I,K)=1.0
           ASYCL(I,K)=1.0
           TAUX(I,K)=TAUCLD(I,K,1)+TAUCLD(I,K,2)
           IF (TAUX(I,K).GT.0.05) THEN
              REFF1=MIN(REFF(I,K,1),130.)
              REFF2=MIN(REFF(I,K,2),20.0)
              W1=(1.-(AIA1+(AIA2+AIA3*REFF1)*REFF1))*TAUCLD(I,K,1)
              W2=(1.-(AWA1+(AWA2+AWA3*REFF2)*REFF2))*TAUCLD(I,K,2)
              SSACL(I,K)=(W1+W2)/TAUX(I,K)
              G1=(AIG1+(AIG2+AIG3*REFF1)*REFF1)*W1
              G2=(AWG1+(AWG2+AWG3*REFF2)*REFF2)*W2
              ASYCL(I,K)=(G1+G2)/(W1+W2)
          ENDIF
         ENDDO
       ENDDO

            DO I=1,MMM
              RR(I,NP+1)=RSIRBM(I)
              RS(I,NP+1)=RSIRBM(I)
              TD(I,NP+1)=0.0
              TT(I,NP+1)=0.0
              TS(I,NP+1)=0.0
            ENDDO

C-----INTEGRATION OVER THE K-DISTRIBUTION FUNCTION

         DO 200 IK=1,NK
c$doacross local(k,i,tausto,tauto,ssato,sssato,asyto,zth,ff,xx,taup,
c$&              sscp,gp,gm1,gm2,gm12p,gm12m,akk,st7,st8,st3,gm3,alf1,
c$&              alf2,all,bll,cll,dll,fll,ell,st2,st4,st1,ONEA,UUU,
c$&              TTT,EMT,UP1,UM1)

          DO 300 K= 1, NP
            DO I= 1, MMM

             TAUSTO=XK(IK)*WH(I,K)+1.0E-8
             TAUTO=TAUSTO+TAUX(I,K)
             SSATO=(SSACL(I,K)*TAUX(I,K))/TAUTO+1.0E-8
             SSATO=MIN(SSATO,0.999999)
             ASYTO=(ASYCL(I,K)*SSACL(I,K)*TAUX(I,K))/(SSATO*TAUTO)

                 zth = 1. / csm(i)

c  delta-eddington scaling of single scattering albedo,
c  optical thickness, and asymmetry factor,
c  K & H eqs(27-29)

                ff  = asyto*asyto
                xx  = 1.-ff*ssato
                taup= tauto*xx
                sscp= ssato*(1.-ff)/xx
                gp  = asyto/(1.+asyto)

c  gamma1, gamma2, and gamma3. see table 2 and eq(26) K & H
c  ssc and gp are the d-s single scattering
c  albedo and asymmetry factor.

                xx  =  3.*gp
                gm1 =  (7. - sscp*(4.+xx))*.25
                gm2 = -(1. - sscp*(4.-xx))*.25
                gm12p=gm1+gm2
                gm12m=gm1-gm2

c  akk is k as defined in eq(25) of K & H

                akk = sqrt(gm12p*gm12m)
                xx  = akk * zth
                st7 = 1. - xx
                st8 = 1. + xx
                st3 = st7 * st8

                if (abs(st3) .lt. 1.E-8) then
                    zth = zth + 0.001
                    xx  = akk * zth
                    st7 = 1. - xx
                    st8 = 1. + xx
                    st3 = st7 * st8
                endif

c  extinction of the direct beam transmission

                td(i,k)  = exp(-taup/zth)

c  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of K & H

                gm3  = (2. - zth*3.*gp)*.25
                alf1 = gm1 - gm3 * gm12m
                alf2 = gm2 + gm3 * gm12m

c  all is last term in eq(21) of K & H
c  bll is last term in eq(22) of K & H

                xx  = akk * 2.
                all = (gm3 - alf2 * zth)*xx*td(i,k)
                bll = (1. - gm3 + alf1*zth)*xx
                xx  = akk * gm3
                cll = (alf2 + xx) * st7
                dll = (alf2 - xx) * st8
                xx  = akk * (1.-gm3)
                fll = (alf1 + xx) * st8
                ell = (alf1 - xx) * st7
                st2 = exp(-akk*taup)
                st4 = st2 * st2
                st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)

c  rr is r-hat of eq(21) of K & H
c  tt is diffuse part of t-hat of eq(22) of K & H

                rr(i,k)=( cll-dll*st4-all*st2)*st1
                tt(i,k)=-((fll-ell*st4)*td(i,k)-bll*st2)*st1
                rr(i,k) = max(rr(i,k),0.)
                tt(i,k) = max(tt(i,k),0.)
                IF(ssato .GT. 0.001) THEN
                   ONEA=1.-ssato
                   XX  = 1.-ssato*asyto
                   UUU = SQRT(XX/ONEA)
                   TTT = SQRT(XX*ONEA*3.)*tauto
                   EMT = EXPMN(-TTT)
                   UP1 = UUU + 1.
                   UM1 = UUU - 1.
                   XX  = UM1*EMT
                   ST1 = 1. / ((UP1+XX) * (UP1-XX))
                   rs(i,k) = UP1*UM1*(1.-EMT*EMT)*ST1
                   ts(i,k) = UUU*4.*EMT*ST1
                ELSE
                   rs(i,k) = 0.0
                   ts(i,k) = EXPMN(-1.66*tauto)
                ENDIF

           ENDDO
 300  CONTINUE

C-----FLUX CALCULATIONS AT EACH LEVEL USING THE TWO-STREAM ADDING METHOD
 
       CALL ADDING (M,RR,TT,TD,RS,TS,FLXDN,FDNDIR,FDNDIF)

       DO K= 1, NP+1
         DO I= 1, MMM
          FLX(I,K) = FLX(I,K)+FLXDN(I,1,K)*HK(IB,IK)
         ENDDO
       ENDDO

        DO I= 1, MMM
          FDIRIR(I) = FDIRIR(I) + FDNDIR(I)*HK(IB,IK)
          FDIFIR(I) = FDIFIR(I) + FDNDIF(I)*HK(IB,IK)
        ENDDO

  200 CONTINUE
  100 CONTINUE
 
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOLUV (M,OH,DP,TAUCLD,REFF,
     *                  CSM,RSUVBM,FLX,FDIRPAR,FDIFPAR)

c      IMPLICIT NONE

      INTEGER NXX,NZ,NADD,NP,NXI,NBB,MM,NBAND
      PARAMETER (NXX=130,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER (NXI=NXX-2)
      PARAMETER (NBB=1,MM=NXI/NBB)
      PARAMETER (NBAND=8)
C************************************************************************
C  COMPUTE SOLAR FLUXES IN THE UV+VISIBLE REGION. THE SPECTRUM IS
C  GROUPED INTO 8 BANDS:
C  
C              BAND     MICROMETER
C
C       UV-C    1.     .175 - .225
C               2.     .225 - .245
C                      .260 - .280
C               3.     .245 - .260
C
C       UV-B    4.     .280 - .295
C               5.     .295 - .310
C               6.     .310 - .320
C      
C       UV-A    7.     .320 - .400
C      
C       PAR     8.     .400 - .700
C
C----- INPUT PARAMETERS:                            UNITS      SIZE
C
C   NUMBER OF SOUNDINGS IN ZONAL DIRECTION (M)       N/D        1
C   NUMBER OF SOUNDINGS IN MERIDIONAL DIRECTION (N)  N/D        1
C   MAXIMUM NUMBER OF SOUNDINGS IN                   N/D        1
C           MERIDIONAL DIRECTION (NDIM)
C   NUMBER OF ATMOSPHERIC LAYERS (NP)                N/D        1
C   LAYER OZONE CONTENT (OH)                      (CM-ATM)STP M*N*NP
C   LAYER PRESSURE THICKNESS (DP)                    MB       M*N*NP
C   CLOUD OPTICAL THICKNESS (TAUCLD)                 N/D      M*NDIM*NP*2
C          INDEX 1 FOR ICE PATICLES
C          INDEX 2 FOR LIQUID PARTICLES
C   EFFECTIVE CLOUD-PARTICLE SIZE (REFF)           MICROMETER M*NDIM*NP*2
C          INDEX 1 FOR ICE PATICLES
C          INDEX 2 FOR LIQUID PARTICLES
C   AEROSOL OPTICAL THICKNESS (TAUAL)                N/D      M*NDIM*NP
C   COSECANT OF THE SOLAR ZENITH ANGLE (CSM)         N/D      M*N
C   UV+PAR SURFACE ALBEDO FOR BEAM                 FRACTION   M*NDIM
C           RADIATION (RSUVBM)
C
C----- OUTPUT (UPDATED) PARAMETERS:
C
C   ALL-SKY NET DOWNWARD FLUX (FLX)                FRACTION   M*NDIM*(NP+1)
C   ALL-SKY DIRECT DOWNWARD PAR FLUX AT
C          THE SURFACE (FDIRPAR)                   FRACTION   M*NDIM
C   ALL-SKY DIFFUSE DOWNWARD PAR FLUX AT
C          THE SURFACE (FDIFPAR)                   FRACTION   M*NDIM
C
C----- NOTE: THE FOLLOWING PARAMETERS MUST BE SPECIFIED BY USERS:
C
C   AEROSOL SINGLE SCATTERING ALBEDO (SSAAL)         N/D        1
C   AEROSOL ASYMMETRY FACTOR (ASYAL)                 N/D        1
C
*
C***********************************************************************
C-----INPUT PARAMETERS

      INTEGER M
      REAL    TAUCLD(M,NP,1),REFF(M,NP,2)
      REAL    OH(M,1),DP(M,1)
      REAL    RSUVBM(1),CSM(1)

C-----OUTPUT (UPDATED) PARAMETER

      REAL    FLX(M,1)
      REAL    FDIRPAR(1),FDIFPAR(1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C-----STATIC PARAMETERS

C      INTEGER NBAND

      REAL    HK(NBAND),XK(NBAND),RY(NBAND)
      REAL    AIG(3),AWG(3)

C-----TEMPORARY ARRAY

      INTEGER I,K,IB,MMM
      REAL    TAURS,TAUOZ,TAUSTO,TAUTO,SSATO,ASYTO
      REAL    TAUX(MM,NP),REFF1,REFF2,G1,G2,ASYCL(MM,NP)
      REAL    TD(MM,NP+1),RR(MM,NP+1),TT(MM,NP+1),
     $        RS(MM,NP+1),TS(MM,NP+1)
      REAL    FLXDN(MM,NP+1),FDNDIR(MM),FDNDIF(MM)

C-----HK IS THE FRACTIONAL EXTRA-TERRESTRIAL SOLAR FLUX.
C     THE SUM OF HK IS 0.47074.

      DATA HK/.00057, .00367, .00083, .00417,
     *        .00600, .00556, .05913, .39081/

C-----XK IS THE OZONE ABSORPTION COEFFICIENT. UNIT: /(CM-ATM)STP

      DATA XK /30.47, 187.2,  301.9,   42.83,
     *         7.09,  1.25,   0.0345,  0.0539/

C-----RY IS THE EXTINCTION COEFFICIENT FOR RAYLEIGH SCATTERING.
C     UNIT: /MB.

      DATA RY /.00604, .00170, .00222, .00132,
     *         .00107, .00091, .00055, .00012/

C-----COEFFICIENTS FOR COMPUTING THE ASYMMETRY FACTOR OF ICE CLOUDS
C     FROM ASYCL=AIG(*,1)+AIG(*,2)*REFF+AIG(*,3)*REFF**2

      DATA AIG/.74625000,.00105410,-.00000264/

C-----COEFFICIENTS FOR COMPUTING THE ASYMMETRY FACTOR OF LIQUID
C     CLOUDS FROM ASYCL=AWG(*,1)+AWG(*,2)*REFF+AWG(*,3)*REFF**2

      DATA AWG/.82562000,.00529000,-.00014866/

      SAVE

C-----INITIALIZE SURFACE REFLECTANCES AND TRANSMITTANCES
CTAO
      MMM=M
CTAO

C-----COMPUTE CLOUD ASYMMETRY FACTOR FOR A MIXTURE OF
C     LIQUID AND ICE PARTICLES.  UNIT OF REFF IS MICROMETERS.


      DO K= 1, NP
        DO I= 1, MMM

           ASYCL(I,K)=1.0

           TAUX(I,K)=TAUCLD(I,K,1)+TAUCLD(I,K,2)
          IF (TAUX(I,K).GT.0.05) THEN

           REFF1=MIN(REFF(I,K,1),130.)
           REFF2=MIN(REFF(I,K,2),20.0)

           G1=(AIG(1)+(AIG(2)+AIG(3)*REFF1)*REFF1)*TAUCLD(I,K,1)
           G2=(AWG(1)+(AWG(2)+AWG(3)*REFF2)*REFF2)*TAUCLD(I,K,2)
           ASYCL(I,K)=(G1+G2)/TAUX(I,K)

          ENDIF
       ENDDO
      ENDDO

         DO I=1,MMM
         RR(I,NP+1)=RSUVBM(I)
         RS(I,NP+1)=RSUVBM(I)
         TD(I,NP+1)=0.0
         TT(I,NP+1)=0.0
         TS(I,NP+1)=0.0
         ENDDO
            
C-----INTEGRATION OVER SPECTRAL BANDS

      DO 100 IB=1,NBAND
c$doacross local(k,i,tausto,tauto,ssato,sssato,asyto,zth,ff,xx,taup,
c$&              sscp,gp,gm1,gm2,gm12p,gm12m,akk,st7,st8,st3,gm3,alf1,
c$&              alf2,all,bll,cll,dll,fll,ell,st2,st4,st1,ONEA,UUU,
c$&              TTT,EMT,UP1,UM1,taurs,tauoz)
       DO 300 K= 1, NP
         DO I= 1, MMM

C-----COMPUTE OZONE AND RAYLEIGH OPTICAL THICKNESSES

          TAURS=RY(IB)*DP(I,K)
          TAUOZ=XK(IB)*OH(I,K)
 
C-----COMPUTE TOTAL OPTICAL THICKNESS, SINGLE SCATTERING ALBEDO,
C     AND ASYMMETRY FACTOR

          TAUSTO=TAURS+TAUOZ+1.0E-8

C-----COMPUTE REFLECTANCE AND TRANSMITTANCE

           TAUTO=TAUSTO+TAUX(I,K)
           SSATO=(TAURS+TAUX(I,K))/TAUTO+1.0E-8
           SSATO=MIN(SSATO,0.999999)
           ASYTO=(ASYCL(I,K)*TAUX(I,K))/(SSATO*TAUTO)
                 zth = 1. / csm(i)

c  delta-eddington scaling of single scattering albedo,
c  optical thickness, and asymmetry factor,
c  K & H eqs(27-29)

                ff  = asyto*asyto
                xx  = 1.-ff*ssato
                taup= tauto*xx
                sscp= ssato*(1.-ff)/xx
                gp  = asyto/(1.+asyto)

c  gamma1, gamma2, and gamma3. see table 2 and eq(26) K & H
c  ssc and gp are the d-s single scattering
c  albedo and asymmetry factor.

                xx  =  3.*gp
                gm1 =  (7. - sscp*(4.+xx))*.25
                gm2 = -(1.   - sscp*(4.-xx))*.25
                gm12p=gm1+gm2
                gm12m=gm1-gm2

c  akk is k as defined in eq(25) of K & H

                akk = sqrt(gm12p*gm12m)

                xx  = akk * zth
                st7 = 1. - xx
                st8 = 1. + xx
                st3 = st7 * st8

                if (abs(st3) .lt. 1.E-8) then
                    zth = zth + 0.001
                    xx  = akk * zth
                    st7 = 1. - xx
                    st8 = 1. + xx
                    st3 = st7 * st8
                endif

c  extinction of the direct beam transmission

                td(i,k)  = exp(-taup/zth)

c  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of K & H

                gm3  = (2. - zth*3.*gp)*.25
                alf1 = gm1 - gm3 * gm12m
                alf2 = gm2 + gm3 * gm12m

c  all is last term in eq(21) of K & H
c  bll is last term in eq(22) of K & H

                xx  = akk * 2.
                all = (gm3 - alf2 * zth    )*xx*td(i,k)
                bll = (1. - gm3 + alf1*zth)*xx

                xx  = akk * gm3
                cll = (alf2 + xx) * st7
                dll = (alf2 - xx) * st8

                xx  = akk * (1.-gm3)
                fll = (alf1 + xx) * st8
                ell = (alf1 - xx) * st7
                st2 = exp(-akk*taup)
                st4 = st2 * st2
                st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)

c  rr is r-hat of eq(21) of K & H
c  tt is diffuse part of t-hat of eq(22) of K & H

                rr(i,k)=( cll-dll*st4-all*st2)*st1
                tt(i,k)=-((fll-ell*st4)*td(i,k)-bll*st2)*st1
                rr(i,k) = max(rr(i,k),0.)
                tt(i,k) = max(tt(i,k),0.)
                IF(ssato .GT. 0.001) THEN
                   ONEA=1.-ssato
                   XX  = 1.-ssato*asyto
                   UUU = SQRT(XX/ONEA)
                   TTT = SQRT(XX*ONEA*3.)*tauto
                   EMT = EXPMN(-TTT)
                   UP1 = UUU + 1.
                   UM1 = UUU - 1.
                   XX  = UM1*EMT
                   ST1 = 1. / ((UP1+XX) * (UP1-XX))
                   rs(i,k) = UP1*UM1*(1.-EMT*EMT)*ST1
                   ts(i,k) = UUU*4.*EMT*ST1
                ELSE
                   rs(i,k) = 0.0
                   ts(i,k) = EXPMN(-1.66*tauto)
                ENDIF

        ENDDO
 300  CONTINUE

C-----FLUX CALCULATIONS

       CALL ADDING (MMM,RR,TT,TD,RS,TS,FLXDN,FDNDIR,FDNDIF)

       DO K= 1, NP+1
         DO I= 1, MMM
          FLX(I,K)=FLX(I,K)+FLXDN(I,K)*HK(IB)
         ENDDO
       ENDDO

       IF(IB.EQ.8) THEN
          DO I=1,MMM
           FDIRPAR(I)=FDNDIR(I)*HK(IB)
           FDIFPAR(I)=FDNDIF(I)*HK(IB)
         ENDDO
       ENDIF

 100  CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION EXPMN(FIN)

c      IMPLICIT NONE
C*******************************************************************
C COMPUTE EXPONENTIAL FOR ARGUMENTS IN THE RANGE 0> FIN > -10.

      REAL    ONE,EXPMIN,E1,E2,E3,E4
      PARAMETER (ONE=1.0, EXPMIN=-10.0)
      PARAMETER (E1=1.0,        E2=-2.507213E-1)
      PARAMETER (E3=2.92732E-2, E4=-3.827800E-3)

      REAL    FIN,EXPMN

      IF (FIN .LT. EXPMIN) FIN = EXPMIN
      EXPMN = ((E4*FIN + E3)*FIN+E2)*FIN+E1
      EXPMN = EXPMN * EXPMN
      EXPMN = ONE / (EXPMN * EXPMN)

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ADDING (M,RR,TT,TD,RS,TS,FLXDN,FDNDIR,FDNDIF)

      INTEGER NX,NZ,NADD,NP,NXI,NBB,MM
      PARAMETER (NX=130,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER (NXI=NX-2)
      PARAMETER (NBB=1,MM=NXI/NBB)

C*********************************************************************
C  COMPUTE UPWARD AND DOWNWARD FLUXES USING A TWO-STREAM ADDING METHOD
C  COMPUTATIONS FOLLOW EQUATIONS (3)-(5) OF CHOU (1992, JAS)
C  INPUT PARAMETERS:
C     NP: TOTAL NUMBER OF LAYERS 
C     RR:  REFLECTION OF A LAYER ILLUMINATED BY BEAM RADIATION
C     TT:  DIFFUSE TRANSMISSION OF A LAYER ILLUMINATED BY BEAM RADIATION
C     TD: DIRECT BEAM TRANMSSION
C     TS: TRANSMISSION OF A LAYER ILLUMINATED BY DIFFUSE RADIATION
C     RS: REFLECTION OF A LAYER ILLUMINATED BY DIFFUSE RADIATION
C
C  OUTPUT PARAMETERS:
C     FLXDN:  NET DOWNWARD FLUXES
C     FDNDIR: SURFACE DIRECT DOWNWARD FLUX
C     FDNDIF: SURFACE DIFFUSE DOWNWARD FLUX
C*********************************************************************

C-----INPUT PARAMETERS

      INTEGER M
      REAL  RR(M,1),TT(M,1),TD(M,1),RS(M,1),
     $        TS(M,1)

C-----OUTPUT PARAMETERS

      REAL    FLXDN(M,1),FDNDIR(1),FDNDIF(1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C-----TEMPORARY ARRAY

      INTEGER I,K,MMM
      REAL    DENM,XX,FUPDIF
      REAL    RSSAB(MM,NP+1),RABX(MM,NP+1),RSABX(MM,NP+1)
      REAL    TBAB(MM,NP+1),TAB(MM,NP+1)

      SAVE

C-----LAYERS ARE ADDED ONE AT A TIME, GOING DOWN FROM THE TOP LAYER,
C     TBAB IS THE COMPOSITE TRANSMITTANCE ILLUMINATED BY BEAM RADIATION
C     TAB IS THE COMPOSITE DIFFUSE TRANSMITTANCE ILLUMINATED BY
C         BEAM RADIATION
C     RSSAB IS THE COMPOSITE REFLECTANCE ILLUMINATED FROM BELOW
C         BY DIFFUSE RADIATION
C     TAB AND RSSAB ARE COMPUTED FROM EQS. (4B) AND (3B) OF CHOU
CTAO
      MMM=M
CTAO
        DO I= 1, MMM
          TBAB(I,1)  = TD(I,1)
          TAB(I,1)   = TT(I,1)
          RSSAB(I,1) = RS(I,1)
        ENDDO

      DO K= 2, NP
       KM=K-1
        DO I= 1, MMM
                 DENM = TS(I,K)/( 1.-RSSAB(I,KM)*RS(I,K) )
          TBAB(I,K) = TBAB(I,KM)*TD(I,K)
          TAB(I,K)  = TBAB(I,KM)*TT(I,K)
     1                + (TBAB(I,KM)*RSSAB(I,KM)
     2                * RR(I,K)+TAB(I,KM))*DENM
          RSSAB(I,K)= RS(I,K)+TS(I,K)*RSSAB(I,KM)*DENM
        ENDDO
      ENDDO

C-----LAYERS ARE ADDED ONE AT A TIME, GOING UP
C     RABX IS THE COMPOSITE REFLECTANCE ILLUMINATED BY BEAM RADIATION
C     RSABX IS THE COMPOSITE REFLECTANCE ILLUMINATED FROM ABOVE
C         BY DIFFUSE RADIATION
C     RABX AND RSABX ARE COMPUTED FROM EQS. (4A) AND (3A) OF CHOU
 
       DO I= 1, MMM
         RABX(I,NP+1) = RR(I,NP+1)
         RSABX(I,NP+1)= RS(I,NP+1)
       ENDDO

      DO K= NP, 1, -1
        KP=K+1
        DO I= 1, MMM
                 DENM  = TS(I,K)/( 1.-RS(I,K)*RSABX(I,KP) )
          RABX(I,K)  = RR(I,K)+(TD(I,K)*RABX(I,KP)
     *                 + TT(I,K)*RSABX(I,KP))*DENM
          RSABX(I,K) = RS(I,K)+TS(I,K)*RSABX(I,KP)*DENM
        ENDDO
      ENDDO
 
C-----COMPUTE FLUXES FOLLOWING EQ (5) OF CHOU (1992)
 
C     FDNDIR IS THE DIRECT  DOWNWARD FLUX
C     FDNDIF IS THE DIFFUSE DOWNWARD FLUX
C     FUPDIF IS THE DIFFUSE UPWARD FLUX

      DO K=2,NP+1
        DO I=1, MMM
                DENM  = 1./(1.- RSSAB(I,K-1)*RSABX(I,K))
         FDNDIR(I)  = TBAB(I,K-1)
                  XX  = TBAB(I,K-1)*RABX(I,K)
         FDNDIF(I)  = (XX*RSSAB(I,K-1)+TAB(I,K-1))*DENM
         FUPDIF  = (XX+TAB(I,K-1)*RSABX(I,K))*DENM
         FLXDN(I,K) = FDNDIR(I)+FDNDIF(I)-FUPDIF
        ENDDO
      ENDDO

        DO I=1, MMM
         FLXDN(I,1) = 1.0-RABX(I,1)
        ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLXCO2(M,NP,SWC,SWH,CSM,DF)

C*****************************************************************

C-----COMPUTE THE REDUCTION OF CLEAR-SKY DOWNWARD SOLAR FLUX
C     DUE TO CO2 ABSORPTION.

c      IMPLICIT NONE

C-----INPUT PARAMETERS


      INTEGER M,NP
      REAL    CSM(1),SWC(M,NP+1),SWH(M,NP+1),CAH(22,19)

C-----OUTPUT (UNDATED) PARAMETER

      REAL    DF(M,NP+1)

      INTEGER IRADAVE
      COMMON/IPTIONR/ IRADAVE

C-----TEMPORARY ARRAY

      INTEGER I,J,K,IC,IW,MMM
      REAL    XX,CLOG,WLOG,DC,DW,X1,X2,Y2

C********************************************************************
C-----INCLUDE CO2 LOOK-UP TABLE

      include "cah.dat"		!cray

      SAVE
 
C********************************************************************
C-----TABLE LOOK-UP FOR THE REDUCTION OF CLEAR-SKY SOLAR
C     RADIATION DUE TO CO2. THE FRACTION 0.0343 IS THE
C     EXTRATERRESTRIAL SOLAR FLUX IN THE CO2 BANDS.
CTAO
      MMM=M
CTAO
      XX=1./.3
      np1=np+1
      DO K= 2, NP1
        DO I= 1, MMM
          CLOG=LOG10(SWC(I,K)*CSM(I))
          WLOG=LOG10(SWH(I,K)*CSM(I))
          IC=INT( (CLOG+3.15)*XX+1.)
          IW=INT( (WLOG+4.15)*XX+1.)
          IC=MAX(2,IC)
          IF(IW.LT.2)IW=2
          IF(IC.GT.22)IC=22
          IF(IW.GT.19)IW=19     
          DC=CLOG-FLOAT(IC-2)*.3+3.
          DW=WLOG-FLOAT(IW-2)*.3+4.   
          X1=CAH(1,IW-1)+(CAH(1,IW)-CAH(1,IW-1))*XX*DW
          X2=CAH(IC-1,IW-1)+(CAH(IC-1,IW)-CAH(IC-1,IW-1))*XX*DW
          Y2=X2+(CAH(IC,IW-1)-CAH(IC-1,IW-1))*XX*DC
          IF (X1.LT.Y2) X1=Y2
          DF(I,K)=DF(I,K)+0.0343*(X1-Y2)
        ENDDO     
      ENDDO      

      RETURN
      END


C********* U1/FRMDC/RADMOD/CLOUD.F *****************                          
C-----LR,LT AND LR2 ARE FOR AN ASYMMETRY FACTOR 0F .843
      BLOCK DATA                                                                
      COMMON/SCO2/CAH(22,19)
      COMMON/O3/ O1(26,21),O2(26,21)
      INTEGER I,IP,IW,J
C                                                                               
      DATA ((O1(IP,IW),IW=1,21),IP=1,3)/                                        
     &-0.3000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,-0.5000E-05,              
     &-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.4000E-05,-0.3000E-05,              
     &-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05, 0.3000E-05,              
     & 0.3000E-05, 0.2000E-05, 0.1000E-05, 0.0000E+00,-0.1000E-05,              
     &-0.3000E-05,-0.3000E-05,-0.3000E-05,-0.3000E-05,-0.4000E-05,              
     &-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.4000E-05,              
     &-0.3000E-05,-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05,              
     & 0.3000E-05, 0.3000E-05, 0.2000E-05, 0.1000E-05, 0.0000E+00,              
     &-0.1000E-05,-0.3000E-05,-0.3000E-05,-0.3000E-05,-0.4000E-05,              
     &-0.4000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,              
     &-0.4000E-05,-0.3000E-05,-0.1000E-05, 0.0000E+00, 0.1000E-05,              
     & 0.2000E-05, 0.3000E-05, 0.3000E-05, 0.2000E-05, 0.1000E-05,              
     & 0.0000E+00,-0.1000E-05,-0.3000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=4,6)/                                        
     &-0.3000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,-0.5000E-05,              
     &-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.4000E-05,-0.3000E-05,              
     &-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05, 0.3000E-05,              
     & 0.3000E-05, 0.2000E-05, 0.1000E-05, 0.0000E+00,-0.1000E-05,              
     &-0.3000E-05,-0.3000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,              
     &-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.4000E-05,              
     &-0.3000E-05,-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05,              
     & 0.3000E-05, 0.3000E-05, 0.2000E-05, 0.1000E-05, 0.0000E+00,              
     &-0.1000E-05,-0.3000E-05,-0.2000E-05,-0.3000E-05,-0.4000E-05,              
     &-0.4000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,              
     &-0.4000E-05,-0.3000E-05,-0.1000E-05, 0.0000E+00, 0.1000E-05,              
     & 0.2000E-05, 0.3000E-05, 0.3000E-05, 0.3000E-05, 0.1000E-05,              
     & 0.0000E+00,-0.1000E-05,-0.2000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=7,9)/                                        
     &-0.2000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,-0.5000E-05,              
     &-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.4000E-05,-0.3000E-05,              
     &-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05, 0.3000E-05,              
     & 0.3000E-05, 0.3000E-05, 0.2000E-05, 0.0000E+00,-0.1000E-05,              
     &-0.2000E-05,-0.3000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,              
     &-0.4000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,-0.4000E-05,              
     &-0.3000E-05,-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05,              
     & 0.3000E-05, 0.3000E-05, 0.3000E-05, 0.2000E-05, 0.0000E+00,              
     &-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.3000E-05,-0.3000E-05,              
     &-0.4000E-05,-0.4000E-05,-0.5000E-05,-0.5000E-05,-0.5000E-05,              
     &-0.4000E-05,-0.2000E-05,-0.1000E-05, 0.0000E+00, 0.1000E-05,              
     & 0.2000E-05, 0.3000E-05, 0.3000E-05, 0.3000E-05, 0.2000E-05,              
     & 0.0000E+00,-0.1000E-05,-0.2000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=10,12)/                                      
     &-0.3000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,-0.4000E-05,              
     &-0.4000E-05,-0.5000E-05,-0.4000E-05,-0.3000E-05,-0.2000E-05,              
     &-0.1000E-05, 0.0000E+00, 0.1000E-05, 0.2000E-05, 0.3000E-05,              
     & 0.3000E-05, 0.3000E-05, 0.1000E-05, 0.0000E+00,-0.1000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.3000E-05,-0.3000E-05,-0.3000E-05,              
     &-0.4000E-05,-0.4000E-05,-0.4000E-05,-0.4000E-05,-0.3000E-05,              
     &-0.2000E-05, 0.0000E+00, 0.0000E+00, 0.1000E-05, 0.2000E-05,              
     & 0.3000E-05, 0.3000E-05, 0.2000E-05, 0.1000E-05, 0.0000E+00,              
     &-0.1000E-05,-0.2000E-05,-0.3000E-05,-0.3000E-05,-0.3000E-05,              
     &-0.3000E-05,-0.3000E-05,-0.4000E-05,-0.4000E-05,-0.3000E-05,              
     &-0.3000E-05,-0.2000E-05, 0.0000E+00, 0.0000E+00, 0.1000E-05,              
     & 0.2000E-05, 0.2000E-05, 0.2000E-05, 0.2000E-05, 0.1000E-05,              
     & 0.0000E+00,-0.1000E-05,-0.3000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=13,15)/                                      
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.3000E-05,-0.3000E-05,              
     &-0.3000E-05,-0.3000E-05,-0.3000E-05,-0.2000E-05,-0.1000E-05,              
     & 0.0000E+00, 0.0000E+00, 0.1000E-05, 0.2000E-05, 0.2000E-05,              
     & 0.2000E-05, 0.1000E-05, 0.0000E+00, 0.0000E+00,-0.1000E-05,              
     &-0.3000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.3000E-05,-0.3000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.1000E-05, 0.0000E+00, 0.0000E+00, 0.1000E-05, 0.1000E-05,              
     & 0.1000E-05, 0.1000E-05, 0.1000E-05, 0.0000E+00, 0.0000E+00,              
     &-0.2000E-05,-0.3000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.1000E-05, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.1000E-05, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     &-0.1000E-05,-0.2000E-05,-0.3000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=16,18)/                                      
     &-0.1000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.1000E-05,              
     &-0.1000E-05, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00,-0.1000E-05,-0.2000E-05,              
     &-0.3000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.1000E-05,-0.1000E-05, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,-0.1000E-05,              
     &-0.2000E-05,-0.3000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     &-0.1000E-05,-0.2000E-05,-0.2000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=19,21)/                                      
     &-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.1000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.1000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,              
     & 0.0000E+00, 0.0000E+00,-0.1000E-05,-0.1000E-05,-0.1000E-05,              
     &-0.2000E-05,-0.1000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.1000E-05,              
     &-0.1000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,              
     &-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.1000E-05,              
     &-0.2000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,              
     &-0.1000E-05,-0.2000E-05,-0.2000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=22,24)/                                      
     &-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.2000E-05,              
     &-0.1000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.1000E-05,-0.1000E-05,-0.1000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.1000E-05,-0.2000E-05,-0.1000E-05,-0.1000E-05,              
     &-0.2000E-05,-0.1000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.3000E-05/                                      
      DATA ((O1(IP,IW),IW=1,21),IP=25,26)/                                      
     &-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.2000E-05,              
     &-0.1000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.3000E-05,-0.2000E-05,-0.1000E-05,-0.2000E-05,-0.1000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.1000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,-0.2000E-05,              
     &-0.2000E-05,-0.3000E-05/                                                  
C                                                                               
      DATA ((O2(IP,IW),IW=1,21),IP=1,3)/                                        
     & 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03, 0.9000E-03,              
     & 0.1300E-02, 0.1700E-02, 0.2200E-02, 0.2800E-02, 0.3300E-02,              
     & 0.3700E-02, 0.4200E-02, 0.4500E-02, 0.4900E-02, 0.5100E-02,              
     & 0.5400E-02, 0.5500E-02, 0.5500E-02, 0.5500E-02, 0.5300E-02,              
     & 0.5200E-02, 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.9000E-03, 0.1300E-02, 0.1700E-02, 0.2200E-02, 0.2800E-02,              
     & 0.3300E-02, 0.3700E-02, 0.4200E-02, 0.4500E-02, 0.4900E-02,              
     & 0.5100E-02, 0.5300E-02, 0.5500E-02, 0.5500E-02, 0.5400E-02,              
     & 0.5300E-02, 0.5100E-02, 0.1000E-03, 0.2000E-03, 0.4000E-03,              
     & 0.6000E-03, 0.9000E-03, 0.1300E-02, 0.1700E-02, 0.2200E-02,              
     & 0.2700E-02, 0.3300E-02, 0.3700E-02, 0.4100E-02, 0.4500E-02,              
     & 0.4800E-02, 0.5100E-02, 0.5300E-02, 0.5400E-02, 0.5400E-02,              
     & 0.5400E-02, 0.5200E-02, 0.5100E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=4,6)/                                        
     & 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03, 0.9000E-03,              
     & 0.1300E-02, 0.1700E-02, 0.2200E-02, 0.2700E-02, 0.3200E-02,              
     & 0.3700E-02, 0.4100E-02, 0.4500E-02, 0.4800E-02, 0.5100E-02,              
     & 0.5300E-02, 0.5400E-02, 0.5400E-02, 0.5300E-02, 0.5200E-02,              
     & 0.5000E-02, 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.9000E-03, 0.1300E-02, 0.1700E-02, 0.2200E-02, 0.2700E-02,              
     & 0.3200E-02, 0.3700E-02, 0.4100E-02, 0.4500E-02, 0.4800E-02,              
     & 0.5000E-02, 0.5200E-02, 0.5300E-02, 0.5300E-02, 0.5200E-02,              
     & 0.5000E-02, 0.4800E-02, 0.1000E-03, 0.2000E-03, 0.4000E-03,              
     & 0.6000E-03, 0.9000E-03, 0.1300E-02, 0.1700E-02, 0.2200E-02,              
     & 0.2700E-02, 0.3200E-02, 0.3700E-02, 0.4100E-02, 0.4400E-02,              
     & 0.4700E-02, 0.5000E-02, 0.5100E-02, 0.5200E-02, 0.5200E-02,              
     & 0.5100E-02, 0.4900E-02, 0.4600E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=7,9)/                                        
     & 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03, 0.9000E-03,              
     & 0.1300E-02, 0.1700E-02, 0.2200E-02, 0.2700E-02, 0.3200E-02,              
     & 0.3600E-02, 0.4000E-02, 0.4300E-02, 0.4600E-02, 0.4800E-02,              
     & 0.5000E-02, 0.5000E-02, 0.5000E-02, 0.4800E-02, 0.4600E-02,              
     & 0.4400E-02, 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.9000E-03, 0.1200E-02, 0.1600E-02, 0.2100E-02, 0.2600E-02,              
     & 0.3100E-02, 0.3500E-02, 0.3900E-02, 0.4200E-02, 0.4500E-02,              
     & 0.4700E-02, 0.4800E-02, 0.4800E-02, 0.4700E-02, 0.4600E-02,              
     & 0.4400E-02, 0.4100E-02, 0.1000E-03, 0.2000E-03, 0.4000E-03,              
     & 0.6000E-03, 0.9000E-03, 0.1200E-02, 0.1600E-02, 0.2100E-02,              
     & 0.2500E-02, 0.3000E-02, 0.3400E-02, 0.3800E-02, 0.4100E-02,              
     & 0.4300E-02, 0.4500E-02, 0.4500E-02, 0.4500E-02, 0.4400E-02,              
     & 0.4200E-02, 0.4000E-02, 0.3800E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=10,12)/                                      
     & 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03, 0.8000E-03,              
     & 0.1100E-02, 0.1500E-02, 0.2000E-02, 0.2400E-02, 0.2900E-02,              
     & 0.3200E-02, 0.3600E-02, 0.3800E-02, 0.4000E-02, 0.4100E-02,              
     & 0.4200E-02, 0.4200E-02, 0.4000E-02, 0.3900E-02, 0.3700E-02,              
     & 0.3400E-02, 0.1000E-03, 0.2000E-03, 0.3000E-03, 0.5000E-03,              
     & 0.7000E-03, 0.1100E-02, 0.1400E-02, 0.1800E-02, 0.2300E-02,              
     & 0.2700E-02, 0.3000E-02, 0.3300E-02, 0.3500E-02, 0.3700E-02,              
     & 0.3800E-02, 0.3800E-02, 0.3800E-02, 0.3600E-02, 0.3500E-02,              
     & 0.3300E-02, 0.3100E-02, 0.0000E+00, 0.1000E-03, 0.3000E-03,              
     & 0.4000E-03, 0.7000E-03, 0.9000E-03, 0.1300E-02, 0.1700E-02,              
     & 0.2000E-02, 0.2400E-02, 0.2700E-02, 0.3000E-02, 0.3200E-02,              
     & 0.3300E-02, 0.3400E-02, 0.3400E-02, 0.3300E-02, 0.3300E-02,              
     & 0.3100E-02, 0.3000E-02, 0.2800E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=13,15)/                                      
     & 0.0000E+00, 0.1000E-03, 0.2000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.8000E-03, 0.1100E-02, 0.1400E-02, 0.1800E-02, 0.2100E-02,              
     & 0.2400E-02, 0.2600E-02, 0.2700E-02, 0.2800E-02, 0.2900E-02,              
     & 0.3000E-02, 0.2900E-02, 0.2900E-02, 0.2800E-02, 0.2700E-02,              
     & 0.2500E-02, 0.0000E+00, 0.1000E-03, 0.1000E-03, 0.3000E-03,              
     & 0.4000E-03, 0.6000E-03, 0.9000E-03, 0.1100E-02, 0.1400E-02,              
     & 0.1700E-02, 0.2000E-02, 0.2200E-02, 0.2300E-02, 0.2400E-02,              
     & 0.2500E-02, 0.2600E-02, 0.2600E-02, 0.2600E-02, 0.2500E-02,              
     & 0.2400E-02, 0.2300E-02, 0.0000E+00, 0.0000E+00, 0.1000E-03,              
     & 0.2000E-03, 0.3000E-03, 0.4000E-03, 0.6000E-03, 0.9000E-03,              
     & 0.1100E-02, 0.1400E-02, 0.1600E-02, 0.1800E-02, 0.1900E-02,              
     & 0.2100E-02, 0.2200E-02, 0.2200E-02, 0.2300E-02, 0.2300E-02,              
     & 0.2300E-02, 0.2200E-02, 0.2200E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=16,18)/                                      
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.1000E-03, 0.2000E-03,              
     & 0.3000E-03, 0.4000E-03, 0.6000E-03, 0.8000E-03, 0.1000E-02,              
     & 0.1200E-02, 0.1400E-02, 0.1600E-02, 0.1700E-02, 0.1900E-02,              
     & 0.2000E-02, 0.2100E-02, 0.2100E-02, 0.2100E-02, 0.2100E-02,              
     & 0.2000E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.1000E-03, 0.2000E-03, 0.3000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.7000E-03, 0.9000E-03, 0.1100E-02, 0.1300E-02, 0.1500E-02,              
     & 0.1600E-02, 0.1800E-02, 0.1900E-02, 0.2000E-02, 0.2000E-02,              
     & 0.2000E-02, 0.1900E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.1000E-03, 0.1000E-03, 0.2000E-03,              
     & 0.4000E-03, 0.5000E-03, 0.7000E-03, 0.9000E-03, 0.1100E-02,              
     & 0.1200E-02, 0.1400E-02, 0.1600E-02, 0.1700E-02, 0.1800E-02,              
     & 0.1900E-02, 0.1900E-02, 0.1900E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=19,21)/                                      
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.1000E-03, 0.1000E-03, 0.2000E-03, 0.4000E-03,              
     & 0.5000E-03, 0.7000E-03, 0.9000E-03, 0.1100E-02, 0.1200E-02,              
     & 0.1400E-02, 0.1600E-02, 0.1700E-02, 0.1800E-02, 0.1800E-02,              
     & 0.1800E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.1000E-03, 0.1000E-03,              
     & 0.2000E-03, 0.4000E-03, 0.5000E-03, 0.7000E-03, 0.9000E-03,              
     & 0.1100E-02, 0.1300E-02, 0.1400E-02, 0.1600E-02, 0.1700E-02,              
     & 0.1700E-02, 0.1700E-02, 0.0000E+00,-0.1000E-03, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.1000E-03, 0.2000E-03, 0.3000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.8000E-03, 0.1000E-02, 0.1200E-02, 0.1300E-02, 0.1500E-02,              
     & 0.1600E-02, 0.1700E-02, 0.1700E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=22,24)/                                      
     &-0.1000E-03,-0.1000E-03, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.1000E-03,              
     & 0.2000E-03, 0.3000E-03, 0.5000E-03, 0.7000E-03, 0.9000E-03,              
     & 0.1100E-02, 0.1300E-02, 0.1400E-02, 0.1500E-02, 0.1600E-02,              
     & 0.1600E-02,-0.1000E-03,-0.1000E-03,-0.1000E-03, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.1000E-03, 0.1000E-03, 0.3000E-03, 0.4000E-03, 0.6000E-03,              
     & 0.8000E-03, 0.1000E-02, 0.1200E-02, 0.1400E-02, 0.1500E-02,              
     & 0.1600E-02, 0.1600E-02,-0.1000E-03,-0.1000E-03, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.1000E-03, 0.2000E-03, 0.4000E-03,              
     & 0.5000E-03, 0.7000E-03, 0.9000E-03, 0.1100E-02, 0.1300E-02,              
     & 0.1500E-02, 0.1600E-02, 0.1600E-02/                                      
      DATA ((O2(IP,IW),IW=1,21),IP=25,26)/                                      
     &-0.1000E-03,-0.1000E-03, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.1000E-03, 0.2000E-03, 0.3000E-03, 0.5000E-03, 0.7000E-03,              
     & 0.9000E-03, 0.1100E-02, 0.1300E-02, 0.1500E-02, 0.1600E-02,              
     & 0.1600E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00,              
     & 0.0000E+00, 0.1000E-03, 0.2000E-03, 0.3000E-03, 0.5000E-03,              
     & 0.7000E-03, 0.9000E-03, 0.1100E-02, 0.1300E-02, 0.1500E-02,              
     & 0.1600E-02, 0.1600E-02/                                                  
C                                                                               
      DATA ((CAH(I,J),I=1,22),J= 1, 5)/                                         
     & 0.9923, 0.9922, 0.9921, 0.9920, 0.9916, 0.9910, 0.9899, 0.9882,          
     & 0.9856, 0.9818, 0.9761, 0.9678, 0.9558, 0.9395, 0.9188, 0.8945,          
     & 0.8675, 0.8376, 0.8029, 0.7621, 0.7154, 0.6647, 0.9876, 0.9876,          
     & 0.9875, 0.9873, 0.9870, 0.9864, 0.9854, 0.9837, 0.9811, 0.9773,          
     & 0.9718, 0.9636, 0.9518, 0.9358, 0.9153, 0.8913, 0.8647, 0.8350,          
     & 0.8005, 0.7599, 0.7133, 0.6627, 0.9808, 0.9807, 0.9806, 0.9805,          
     & 0.9802, 0.9796, 0.9786, 0.9769, 0.9744, 0.9707, 0.9653, 0.9573,          
     & 0.9459, 0.9302, 0.9102, 0.8866, 0.8604, 0.8311, 0.7969, 0.7565,          
     & 0.7101, 0.6596, 0.9708, 0.9708, 0.9707, 0.9705, 0.9702, 0.9697,          
     & 0.9687, 0.9671, 0.9647, 0.9612, 0.9560, 0.9483, 0.9372, 0.9221,          
     & 0.9027, 0.8798, 0.8542, 0.8253, 0.7916, 0.7515, 0.7054, 0.6551,          
     & 0.9568, 0.9568, 0.9567, 0.9565, 0.9562, 0.9557, 0.9548, 0.9533,          
     & 0.9510, 0.9477, 0.9428, 0.9355, 0.9250, 0.9106, 0.8921, 0.8700,          
     & 0.8452, 0.8171, 0.7839, 0.7443, 0.6986, 0.6486/                          
      DATA ((CAH(I,J),I=1,22),J= 6,10)/                                         
     & 0.9377, 0.9377, 0.9376, 0.9375, 0.9372, 0.9367, 0.9359, 0.9345,          
     & 0.9324, 0.9294, 0.9248, 0.9181, 0.9083, 0.8948, 0.8774, 0.8565,          
     & 0.8328, 0.8055, 0.7731, 0.7342, 0.6890, 0.6395, 0.9126, 0.9126,          
     & 0.9125, 0.9124, 0.9121, 0.9117, 0.9110, 0.9098, 0.9079, 0.9052,          
     & 0.9012, 0.8951, 0.8862, 0.8739, 0.8579, 0.8385, 0.8161, 0.7900,          
     & 0.7585, 0.7205, 0.6760, 0.6270, 0.8809, 0.8809, 0.8808, 0.8807,          
     & 0.8805, 0.8802, 0.8796, 0.8786, 0.8770, 0.8747, 0.8712, 0.8659,          
     & 0.8582, 0.8473, 0.8329, 0.8153, 0.7945, 0.7697, 0.7394, 0.7024,          
     & 0.6588, 0.6105, 0.8427, 0.8427, 0.8427, 0.8426, 0.8424, 0.8422,          
     & 0.8417, 0.8409, 0.8397, 0.8378, 0.8350, 0.8306, 0.8241, 0.8148,          
     & 0.8023, 0.7866, 0.7676, 0.7444, 0.7154, 0.6796, 0.6370, 0.5897,          
     & 0.7990, 0.7990, 0.7990, 0.7989, 0.7988, 0.7987, 0.7983, 0.7978,          
     & 0.7969, 0.7955, 0.7933, 0.7899, 0.7846, 0.7769, 0.7664, 0.7528,          
     & 0.7357, 0.7141, 0.6866, 0.6520, 0.6108, 0.5646/                          
      DATA ((CAH(I,J),I=1,22),J=11,15)/                                         
     & 0.7515, 0.7515, 0.7515, 0.7515, 0.7514, 0.7513, 0.7511, 0.7507,          
     & 0.7501, 0.7491, 0.7476, 0.7450, 0.7409, 0.7347, 0.7261, 0.7144,          
     & 0.6992, 0.6793, 0.6533, 0.6203, 0.5805, 0.5357, 0.7020, 0.7020,          
     & 0.7020, 0.7019, 0.7019, 0.7018, 0.7017, 0.7015, 0.7011, 0.7005,          
     & 0.6993, 0.6974, 0.6943, 0.6894, 0.6823, 0.6723, 0.6588, 0.6406,          
     & 0.6161, 0.5847, 0.5466, 0.5034, 0.6518, 0.6518, 0.6518, 0.6518,          
     & 0.6518, 0.6517, 0.6517, 0.6515, 0.6513, 0.6508, 0.6500, 0.6485,          
     & 0.6459, 0.6419, 0.6359, 0.6273, 0.6151, 0.5983, 0.5755, 0.5458,          
     & 0.5095, 0.4681, 0.6017, 0.6017, 0.6017, 0.6017, 0.6016, 0.6016,          
     & 0.6016, 0.6015, 0.6013, 0.6009, 0.6002, 0.5989, 0.5967, 0.5932,          
     & 0.5879, 0.5801, 0.5691, 0.5535, 0.5322, 0.5043, 0.4700, 0.4308,          
     & 0.5518, 0.5518, 0.5518, 0.5518, 0.5518, 0.5518, 0.5517, 0.5516,          
     & 0.5514, 0.5511, 0.5505, 0.5493, 0.5473, 0.5441, 0.5393, 0.5322,          
     & 0.5220, 0.5076, 0.4878, 0.4617, 0.4297, 0.3929/                          
      DATA ((CAH(I,J),I=1,22),J=16,19)/                                         
     & 0.5031, 0.5031, 0.5031, 0.5031, 0.5031, 0.5030, 0.5030, 0.5029,          
     & 0.5028, 0.5025, 0.5019, 0.5008, 0.4990, 0.4960, 0.4916, 0.4850,          
     & 0.4757, 0.4624, 0.4441, 0.4201, 0.3904, 0.3564, 0.4565, 0.4565,          
     & 0.4565, 0.4564, 0.4564, 0.4564, 0.4564, 0.4563, 0.4562, 0.4559,          
     & 0.4553, 0.4544, 0.4527, 0.4500, 0.4460, 0.4400, 0.4315, 0.4194,          
     & 0.4028, 0.3809, 0.3538, 0.3227, 0.4122, 0.4122, 0.4122, 0.4122,          
     & 0.4122, 0.4122, 0.4122, 0.4121, 0.4120, 0.4117, 0.4112, 0.4104,          
     & 0.4089, 0.4065, 0.4029, 0.3976, 0.3900, 0.3792, 0.3643, 0.3447,          
     & 0.3203, 0.2923, 0.3696, 0.3696, 0.3696, 0.3696, 0.3696, 0.3696,          
     & 0.3695, 0.3695, 0.3694, 0.3691, 0.3687, 0.3680, 0.3667, 0.3647,          
     & 0.3615, 0.3570, 0.3504, 0.3409, 0.3279, 0.3106, 0.2892, 0.2642/          
C                                                                               
      END 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PTPB
      PARAMETER (NX=130,NY=130,NZ=43)
      PARAMETER (NXI=NX-2,NYJ=NY-2,NZM1=NZ-1)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)

      COMMON/BW1/ PCLTOP(NXI,NYJ),PCLBOT(NXI,NYJ)
C
      DIMENSION TQE(NZ)
      SAVE
      DO 101 I=2,NX-1
        IX=I-1
      DO 100 J=2,NY-1
        JY=J-1
        DO K=2,NZ-1
          TQE(K)=QCL(I,J,K)+QRN(I,J,K)+QCI(I,J,K)+QCS(I,J,K)+QCG(I,J,K)
        END DO
        CALL TTOP (TQE,KTOP,KBASE)
        PCLTOP(IX,JY)=1.E-3*P0(KTOP+1)
           IF ( KTOP.EQ.NZM1) PCLTOP(IX,JY)=1.E-3*P0(KTOP)
        PCLBOT(IX,JY)=1.E-3*P0(KBASE)
           IF (KBASE.EQ. 2) PCLBOT(IX,JY)=1.E-3*P0(KBASE)
         THIN=PCLBOT(IX,JY)-PCLTOP(IX,JY)
           IF(THIN.LT.50.) GO TO 100
  100 CONTINUE
  101 CONTINUE
      RETURN
      END          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TTOP (TWC,KTOP,KBASE)

c      IMPLICIT NONE
      INTEGER NZ,NZP,NZM
      PARAMETER (NZ=43,NZP=NZ+1,NZM=NZ-1)
      INTEGER KBASE,KTOP
      REAL    TWC(1)

C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER K,K1,KTEMP
      SAVE
C
        KTOP=2
        DO 10 K=2,NZM
          IF (TWC(K).LT.1.E-5) GOTO 10
            KTEMP=K
            IF(KTOP.LT.KTEMP) KTOP=KTEMP
   10   CONTINUE
         KBASE=NZM
        DO 20 K1=2,NZM
          K=NZP-K1
          IF (TWC(K).LT.1.E-5) GOTO 20
            KTEMP=K
            IF (KBASE.GT.KTEMP) KBASE=KTEMP
   20   CONTINUE
        IF(KTOP.LT.KBASE)THEN
          KTOP=1
          KBASE=2
        ENDIF
      RETURN
      END   




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FITO3 (NPP1,P,AT,AO)

c      IMPLICIT NONE
C-----------------------------------------------------------------------
C--- INPUTS: NPP1 --- NUMBER OF POINTS IN THE VERTICAL               ---
C---         P    --- PRESSURE AT THE MODEL LEVELS                   ---
C-----------------------------------------------------------------------
C--- OUTPUTS: AT  --- TEMPRETURE AT THE MODEL LEVELS                 ---
C---          AO  --- OZONE AT THE MODEL LEVELS                      ---
C-----------------------------------------------------------------------
      INTEGER N,NPP1
      PARAMETER (N=61)

      REAL    P(1),AT(1),AO(1)

      REAL    PL(n),TA(n),WA(n),OA(n),PA(n)
      REAL    TL(n),WL(n),OL(n)
      REAL    W(n),TAB(3),WK(N,4)
      INTEGER ITAB(3),IOP(2)
C----- IX=        1:TRP; 2:MLS; 3:MLW; 4:SAS; 5:SAW

C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER I,IX,K,LUN,INT
      REAL    Y

      SAVE

      IX=1
      LUN=90+IX
       OPEN( lun,FILE='trp.dat',FORM='FORMATTED',
     1       STATUS='OLD' ) 	!cray
      READ (LUN,51) (PL(I),TL(I),WL(I),OL(I),I=2,N)
   51 FORMAT (0P,2F10.3,1P,2E12.4)
      PL(1)=0.
      TL(1)=TL(2)
      WL(1)=WL(2)
      OL(1)=OL(2)
      DO 25 I=1,N-1
        IP=I+1
        PA(I)=0.5*(PL(I)+PL(IP))
        TA(I)=0.5*(TL(I)+TL(IP))
        WA(I)=0.5*(WL(I)+WL(IP))
        OA(I)=0.5*(OL(I)+OL(IP))
   25 CONTINUE
      PA(N)=PL(N)
      TA(N)=TL(N)
      WA(N)=WL(N)
      OA(N)=OL(N)
      WRITE (6,105)
      WRITE (6,106) (I,PA(I),TA(I),OA(I),WA(I),I=1,N-1)

C-----------------------------------------------------------------------
C--- STARTS INTERPOLATION                                            ---
C-----------------------------------------------------------------------
      IOP(1)=4
      IOP(2)=4
      INT=1
      ITAB(1)=1
      ITAB(2)=0
      ITAB(3)=0
      CALL COEFF (N,PA,TA,W,IOP,INT,WK)
      DO 200 K=1,NPP1-1
        Y=P(K)
       CALL TERP1(N,PA,TA,W,Y,INT,TAB,ITAB)
        AT(K)=TAB(1)
        PRINT *,'K = ',K,'  Y = ',Y,'  AT = ',AT(K)
  200 CONTINUE
      CALL COEFF (N,PA,OA,W,IOP,INT,WK)
      DO 300 K=1,NPP1-1
        Y=P(K)
       CALL TERP1(N,PA,OA,W,Y,INT,TAB,ITAB)
       AO(K)=TAB(1)
       WRITE(6,106) K,Y,AT(K),AO(K)
  300 CONTINUE
  105 FORMAT(5X,'K',7X,'PA',8X,'TA',10X,'O3',8X,'QV')
  106 FORMAT(2X,I4,2F10.3,2E12.4)
      RETURN
      END

      SUBROUTINE TMAXADV (X,SMALL)
C     ****   FIND MAXIMUM VALUES' ROUTINE
      PARAMETER (NX=130,NY=130,NZ=43)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      DIMENSION X(NX,NY,NZ)
      save
CC    ****   SEARCH FOR THE SMALLEST VALUES
      SMALL=0.
c$doacross local(k,j,i,xa), reduction(small)
      DO K=2,KLES
      DO J=2,JLES
      DO I=2,ILES
       XA=X(I,J,K)
       SMALL=MIN(SMALL,XA)
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE TMAX (X,AA,AMINUT,BIG,SMALL,ID)
C     ****   FIND MAXIMUM AND MINIMUM VALUES' ROUTINE
      PARAMETER (NX=130,NY=130,NZ=43)
      CHARACTER*4 AA
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      DIMENSION X(NX,NY,NZ),CONS(5)
      save
      DATA CONS/1.E-2,1.,1.E3,1.E-7,1.E-5/
CC    ****   SEARCH FOR THE LARGEST OR/AND SMALLEST VALUES
      BIG=0.
      SMALL=0.
c$doacross local(k,j,i,xa), reduction(big,small)
      DO K=2,KLES
      DO J=2,JLES
      DO I=2,ILES
       XA=X(I,J,K)
       BIG=MAX(BIG,XA)
       SMALL=MIN(SMALL,XA)
      ENDDO
      ENDDO
      ENDDO
CC    **************************
      DO K=2,KLES
      DO J=2,JLES
      DO I=2,ILES
       IF (X(I,J,K) .GE. BIG) THEN
        ILO1=I-1
        JLO1=J-1
        KLO1=K-1
       ENDIF
       IF (X(I,J,K) .LE. SMALL) THEN
        ILO2=I-1
        JLO2=J-1
        KLO2=K-1
       ENDIF
      ENDDO
      ENDDO
      ENDDO
CC    **************************
      BIG=CONS(ID)*BIG
      SMALL=CONS(ID)*SMALL
      WRITE(6,123) AMINUT,AA,BIG,ILO1,JLO1,KLO1,SMALL,ILO2,JLO2,KLO2
  123 FORMAT(3X,'AMIN=',F5.0,2X,A4,3X,'BIG=',F6.2,2X,'I=',I3,2X,'J=',
     1   I3,2X,'K=',I2,10X,'SMALL=',F6.2,2X,'I=',I3,2X,'J=',I3,2X,'K=',
     2   I2)
      RETURN
      END

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WMICRO (ITAPE)
C     ******   WRITE DATA ON PHYSICAL TAPE   ******
      PARAMETER (NX=130,NY=130,NZ=43,NT=2880,NT2=NT*2,NT5=NT*5)
      PARAMETER (ITT=244)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B2S/ QCS1(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4A/ AK1(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      COMMON/BRH1/ SRO(NZ),QRO(NZ),QC(NZ),QR(NZ),QI(NZ),QS(NZ),QG(NZ),
     1   TQC(NZ),TQR(NZ),TQI(NZ),TQS(NZ),TQG(NZ),QQQ(NT5)
      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      common/dinradn/ qairsfc(nx,ny),pairsfc(nx,ny),thairsf(nx,ny)

      COMMON/SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
      COMMON/SLWAVE/ RSW(NX,NY,NZ),RLW(NX,NY,NZ)
      COMMON/SURFACE/ SUN_4(NX,NY,4)
      save
C     ****************************************
      write(itape) dpt
      write(itape) dqv
      write(itape) qcl
      write(itape) qrn
      write(itape) qci
      write(itape) qcs
      write(itape) qcg
      write(itape) tb
      write(itape) qb
      write(itape) rho
      write(itape) rho1
      write(itape) ta
      write(itape) qa
      write(itape) tqc
      write(itape) tqr
      write(itape) tqi
      write(itape) tqs
      write(itape) tqg


      RETURN
      END

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WDYN (ITAPE)
C     ******   WRITE DATA ON PHYSICAL TAPE   ******
      PARAMETER (NX=130,NY=130,NZ=43,NT=2880,NT2=NT*2,NT5=NT*5)
      PARAMETER (ITT=244)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B2S/ QCS1(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      COMMON/BRH1/ SRO(NZ),QRO(NZ),QC(NZ),QR(NZ),QI(NZ),QS(NZ),QG(NZ),
     1   TQC(NZ),TQR(NZ),TQI(NZ),TQS(NZ),TQG(NZ),QQQ(NT5)
      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      common/dinradn/ qairsfc(nx,ny),pairsfc(nx,ny),thairsf(nx,ny)

      COMMON/SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
      COMMON/SLWAVE/ RSW(NX,NY,NZ),RLW(NX,NY,NZ)
      COMMON/SURFACE/ SUN_4(NX,NY,4)
      save
C     ****************************************
      write(itape) ak
      write(itape) u
      write(itape) v
      write(itape) w
      write(itape) rsw
      write(itape) rlw
      write(itape) ub
      write(itape) vb
      write(itape) pi
      write(itape) p0
      write(itape) wb
      write(itape) wbt
      write(itape) fd
      write(itape) fe
      write(itape) aq1t
      write(itape) aq2t
      write(itape) aq1zt
      write(itape) aq2zt
      write(itape) suw
      write(itape) svw
      write(itape) swt
      write(itape) swq
      write(itape) sun_4
      write(itape) tairsfc
      write(itape) qairsfc
      write(itape) pairsfc
      write(itape) thairsf

      RETURN
      END

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WSTAT (ITAPE)
C     ******   WRITE DATA ON PHYSICAL TAPE   ******
      PARAMETER (NX=130,NY=130,NZ=43,NT=2880,NT2=2*NT)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      COMMON/B8/ SMF(NZ),SMU(NZ),SMD(NZ),STF(NZ),STU(NZ),STD(NZ),
     1   SQF(NZ),SQU(NZ),SQD(NZ),SDT(NZ),SDQ(NZ),CLD(NZ),STV(NZ)
      COMMON/B9/ SUT1(NZ),SUC1(NZ),SUN1(NZ),SUU1(NZ),SUD1(NZ),AUB(NZ),
     1   SSU(NZ),SUT2(NZ),SUC2(NZ),SUN2(NZ),SUU2(NZ),SUD2(NZ),AUC1(NZ),
     2   AUN1(NZ),AUU1(NZ),AUD1(NZ),AUC2(NZ),AUN2(NZ),AUU2(NZ),AUD2(NZ),
     3   SVT1(NZ),SVC1(NZ),SVN1(NZ),SVU1(NZ),SVD1(NZ),AVB(NZ),
     4   SSV(NZ),SVT2(NZ),SVC2(NZ),SVN2(NZ),SVU2(NZ),SVD2(NZ),AVC1(NZ),
     5   AVN1(NZ),AVU1(NZ),AVD1(NZ),AVC2(NZ),AVN2(NZ),AVU2(NZ),AVD2(NZ),
     6   STF1(NZ),STU1(NZ),STD1(NZ),SQF1(NZ),SQU1(NZ),SQD1(NZ),
     7   STF2(NZ),STU2(NZ),STD2(NZ),SQF2(NZ),SQU2(NZ),SQD2(NZ),CLDU(NZ),
     8   CLDD(NZ),CLDN(NZ),SMN(NZ),STN1(NZ),STN2(NZ),STN(NZ),SQN1(NZ),
     9   SQN2(NZ),SQN(NZ)
      common/bb6/ tls(nz),qls(nz),tls1(nz),tls2(nz),qls1(nz),qls2(nz),
     1   tls3(nz),tls4(nz),qls3(nz),qls4(nz),sft(nz),sfq(nz)
      COMMON/BCS/ S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     1 S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     2 S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     3 SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     4 SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     5 SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCSS/ SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     1 SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     2 SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     3 SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     4 SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     5 SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)
      COMMON/BSFC/ tsfc_1(nx,ny), qsfc_1(nx,ny)
      COMMON/BSTS/ THOM(NZ,4,7),TDW(NZ,4,7),TMLT(NZ,4,7),SAUT(NZ,4,7),
     1 SACI(NZ,4,7),SACW(NZ,4,7),RACI(NZ,4,7),TACR(NZ,4,7),RAUT(NZ,4,7),
     2 RACW(NZ,4,7),SFW(NZ,4,7),SFI(NZ,4,7),GACS(NZ,4,7),GACW(NZ,4,7),
     3 GACI(NZ,4,7),GACR(NZ,4,7),GWET(NZ,4,7),GAUT(NZ,4,7),RACS(NZ,4,7),
     4 SACR(NZ,4,7),GFR(NZ,4,7),SMLT(NZ,4,7),GMLT(NZ,4,7),SDEP(NZ,4,7),
     5 SSUB(NZ,4,7),GSUB(NZ,4,7),PERN(NZ,4,7),D3RI(NZ,4,7),D3IR(NZ,4,7),
     6 D2SR(NZ,4,7),D2RS(NZ,4,7),GDRY(NZ,4,7),COC(NZ,4,7),COE(NZ,4,7),
     7 SMF0(NZ,4,7),QC0(NZ,4,7),QR0(NZ,4,7),QI0(NZ,4,7),QS0(NZ,4,7),
     8 QG0(NZ,4,7),SQC0(NZ,4,7),SQR0(NZ,4,7),SQI0(NZ,4,7),SQS0(NZ,4,7),
     9 SQG0(NZ,4,7),ERNS(NZ,4,7),WGRS(NZ,4,7),QSWS(NZ,4,7),TB0(NZ,4),
     1 QB0(NZ,4)
      COMMON/BSTS1/ TUT1(NZ,4,7),TUT2(NZ,4,7),TVT1(NZ,4,7),TVT2(NZ,4,7),
     1 TSTF(NZ,4,7),TSTF1(NZ,4,7),TSTF2(NZ,4,7),TSQF(NZ,4,7),
     2 QQQ(NZ,4,7),TSQF1(NZ,4,7),TSQF2(NZ,4,7),TSQQ(NZ,4,7),
     3 TSQQ1(NZ,4,7)
      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNTL(NZ,4,7),
     1 SNQH(NZ,4,7),SNQV(NZ,4,7),SNQD(NZ,4,7),SNQL(NZ,4,7),
     2 SNTL1(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7),
     3 SNQL1(NZ,4,7)
      COMMON/BSTS3/ QV0(NZ,4,7),TT0(NZ,4,7),SQV0(NZ,4,7),STT0(NZ,4,7),
     1 SGPT(NZ,4,7),SSGPT(NZ,4,7),SNQHD(NZ,4,7),SNQVD(NZ,4,7),
     2 Q1T(NZ,4,7),SNHDH(NZ,4,7),SQHDT(NZ,4,7),SQVDT(NZ,4,7)
      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7),SQTDT(NZ,4,7),SQHL(NZ,4,7)
      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NX,NY),ICS5(NX,NY,4),
     1  IBZ(NX,NY,4)
      COMMON/BSTS20/ OTHERT_ADD(NZ,4,7),OTHERQ_ADD(NZ,4,7)
      common/bsts40/ fcld(nz,4,7)
      COMMON/BTT/ S1(16,NZ,14),SN1(5,NZ,14)
      common/gbs/ tlsw(nz),qlsw(nz),ttlsw(nz),tqlsw(nz),smfff(nz),
     1   smffff(nz),tcof(nz),textra(nz),qextra(nz)
      save
C     ****************************************
      write(itape) s1
      write(itape) sn1
      write(itape) s9
      write(itape) sn9
      write(itape) s10
      write(itape) sn10
      write(itape) s11
      write(itape) sn11
      write(itape) s12
      write(itape) sn12
      write(itape) s13
      write(itape) sn13
      write(itape) s14
      write(itape) sn14
      write(itape) s15
      write(itape) sn15
      write(itape) s16
      write(itape) sn16
      write(itape) s17
      write(itape) sn17
      write(itape) s18
      write(itape) sn18
      write(itape) s19
      write(itape) sn19
      write(itape) s20
      write(itape) sn20
      write(itape) s21
      write(itape) sn21
      write(itape) scnt
      write(itape) sncnt
      write(itape) ss9
      write(itape) ssn9
      write(itape) ss10
      write(itape) ssn10
      write(itape) ss11
      write(itape) ssn11
      write(itape) ss12
      write(itape) ssn12
      write(itape) ss13
      write(itape) ssn13
      write(itape) ss14
      write(itape) ssn14
      write(itape) ss15
      write(itape) ssn15
      write(itape) ss16
      write(itape) ssn16
      write(itape) ss17
      write(itape) ssn17
      write(itape) ss18
      write(itape) ssn18
      write(itape) ss19
      write(itape) ssn19
      write(itape) ss20
      write(itape) ssn20
      write(itape) ss21
      write(itape) ssn21
      write(itape) sscnt
      write(itape) ssncnt
      write(itape) tls
      write(itape) qls
      write(itape) cld
      write(itape) cldu
      write(itape) cldd
      write(itape) smf
      write(itape) smu
      write(itape) smd
      write(itape) smn
      write(itape) sft
      write(itape) sfq
      write(itape) stf1
      write(itape) stu1
      write(itape) std1
      write(itape) stn1
      write(itape) stf2
      write(itape) stu2
      write(itape) std2
      write(itape) stn2
      write(itape) stf
      write(itape) stu
      write(itape) std
      write(itape) stn
      write(itape) sqf1
      write(itape) squ1
      write(itape) sqd1
      write(itape) sqn1
      write(itape) sqf2
      write(itape) squ2
      write(itape) sqd2
      write(itape) sqn2
      write(itape) sqf
      write(itape) squ
      write(itape) sqd
      write(itape) sqn
      write(itape) sc
      write(itape) se
      write(itape) scu1
      write(itape) sed1
      write(itape) aub
      write(itape) avb
      write(itape) sut1
      write(itape) svt1
      write(itape) suu1
      write(itape) sud1
      write(itape) svu1
      write(itape) svd1
      write(itape) suu2
      write(itape) sud2
      write(itape) svu2
      write(itape) svd2
      write(itape) auu1
      write(itape) aud1
      write(itape) avu1
      write(itape) avd1
      write(itape) auu2
      write(itape) aud2
      write(itape) avu2
      write(itape) avd2
      write(itape) thom
      write(itape) tdw
      write(itape) tmlt
      write(itape) saut
      write(itape) saci
      write(itape) sacw
      write(itape) raci
      write(itape) tacr
      write(itape) raut
      write(itape) racw
      write(itape) sfw
      write(itape) sfi
      write(itape) gacs
      write(itape) gacw
      write(itape) gaci
      write(itape) gacr
      write(itape) gwet
      write(itape) gaut
      write(itape) racs
      write(itape) sacr
      write(itape) gfr
      write(itape) smlt
      write(itape) gmlt
      write(itape) sdep
      write(itape) ssub
      write(itape) gsub
      write(itape) pern
      write(itape) d3ri
      write(itape) d3ir
      write(itape) d2sr
      write(itape) d2rs
      write(itape) gdry
      write(itape) coc
      write(itape) coe
      write(itape) smf0
      write(itape) qc0
      write(itape) qr0
      write(itape) qi0
      write(itape) qs0
      write(itape) qg0
      write(itape) sqc0
      write(itape) sqr0
      write(itape) sqi0
      write(itape) sqs0
      write(itape) sqg0
      write(itape) erns
      write(itape) wgrs
      write(itape) qsws
      write(itape) srsw
      write(itape) srlw
      write(itape) tut1
      write(itape) tut2
      write(itape) tvt1
      write(itape) tvt2
      write(itape) tstf
      write(itape) tstf1
      write(itape) tstf2
      write(itape) tsqf
      write(itape) tsqf1
      write(itape) tsqf2
      write(itape) tsqq
      write(itape) tsqq1
      write(itape) fcld
      write(itape) othert_add
      write(itape) otherq_add
      write(itape) snth
      write(itape) sntv
      write(itape) sntd
      write(itape) snqh
      write(itape) snqv
      write(itape) snqd
      write(itape) snhh
      write(itape) snhv
      write(itape) snhd
      write(itape) tb0
      write(itape) qb0
      write(itape) tsfc_1
      write(itape) qsfc_1
      write(itape) tlsw
      write(itape) qlsw
      write(itape) ttlsw
      write(itape) tqlsw
      write(itape) tcof
      write(itape) textra
      write(itape) qextra

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WTAPRI(ifile)
      PARAMETER (NX=130,NY=130,NT=2880,NT2=NT*2)
      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      WRITE(ifile) RI
      RETURN
      END

C*******************************************************************************

      SUBROUTINE RESTART (FLAG)
C     ******   WRITE DATA ON PHYSICAL TAPE   ******
      PARAMETER (NX=130,NY=130,NZ=43,NT=2880,ITT=244)
      CHARACTER*(*) FLAG
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B1TQ/ DPT(NX,NY,NZ),DQV(NX,NY,NZ)
      COMMON/B2CR/ QCL1(NX,NY,NZ),QRN1(NX,NY,NZ)
      COMMON/B2IG/ QCI1(NX,NY,NZ),QCG1(NX,NY,NZ)
      COMMON/B2S/ QCS1(NX,NY,NZ)
      COMMON/B2TQ/ DPT1(NX,NY,NZ),DQV1(NX,NY,NZ)
      COMMON/B3A/ AK(NX,NY,NZ)
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B3WP/ W(NX,NY,NZ)
      COMMON/B4A/ AK1(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),Z1(NZ),Z2(NZ),Z3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
      COMMON/BRH1/ SRRO(NZ),QRRO(NZ),SQC(NZ),SQR(NZ),SQI(NZ),SQS(NZ),
     1   SQG(NZ),STQC(NZ),STQR(NZ),STQI(NZ),STQS(NZ),STQG(NZ),
     2   TQC(NT),TQR(NT),TQI(NT),TQS(NT),TQG(NT)
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1   NRAN,KT1,KT2
      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      common/dinradn/ qairsfc(nx,ny),pairsfc(nx,ny),thairsf(nx,ny)
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
      COMMON/ISFCRI/ IRICONT
      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/sfcri/ ri180(nx,ny),riold180(nx,ny)
      common/sfcten/tsdt,qsdt,thsdt,psdt
      COMMON/SLTM/ RLAT,RMONTH,RIDAY,HRL,SO0,COSZ,ICOSZ,TERMAN,RSFC
      COMMON/SLWAVE/ RSW(NX,NY,NZ),RLW(NX,NY,NZ)
      COMMON/SUN/ ISUN
      COMMON/SFLUXS/ SUW(NX,NY),SVW(NX,NY),SWT(NX,NY),SWQ(NX,NY)
      save
C     ****************************************

      IF (FLAG .EQ. 'READ') THEN
      
         READ (7) DPT
         READ (7) DQV
         READ (7) QCL
         READ (7) QRN
         READ (7) QCI
         READ (7) QCS
         READ (7) QCG
         READ (7) U
         READ (7) V
         READ (7) W
         READ (7) AK
         READ (7) DPT1
         READ (7) DQV1
         READ (7) QCL1
         READ (7) QRN1
         READ (7) QCI1
         READ (7) QCS1
         READ (7) QCG1
         READ (7) UU1
         READ (7) VV1
         READ (7) WW1
         READ (7) AK1
         READ (7) UMD
         READ (7) VMD
         READ (7) WMD
         READ (7) RSW
         READ (7) RLW
         READ (7) UB
         READ (7) VB
         READ (7) UB1
         READ (7) VB1
         READ (7) TB
         READ (7) QB
         READ (7) RHO
         READ (7) RHO1
         READ (7) TA
         READ (7) QA
         READ (7) TA1
         READ (7) QA1
         READ (7) PI
         READ (7) P0
         READ (7) WB
         READ (7) WBT
         READ (7) UBT
         READ (7) VBT
         READ (7) FD
         READ (7) FE
         READ (7) SQC
         READ (7) SQR
         READ (7) SQI
         READ (7) SQS
         READ (7) SQG
         READ (7) STQC
         READ (7) STQR
         READ (7) STQI
         READ (7) STQS
         READ (7) STQG
         READ (7) Q1T
         READ (7) Q2T
         READ (7) AQ1T
         READ (7) AQ2T
         READ (7) AQ1ZT
         READ (7) AQ2ZT
         READ (7) TAIRSFC
         READ (7) QAIRSFC
         READ (7) THAIRSF
         READ (7) PAIRSFC
         READ (7) TSDT
         READ (7) QSDT
         READ (7) THSDT
         READ (7) PSDT
         READ (7) ISUN
         READ (7) RLAT
         READ (7) RMONTH
         READ (7) RIDAY
         READ (7) IRICONT
         READ (7) HRL
         READ (7) SO0
         READ (7) COSZ
         READ (7) ICOSZ
         READ (7) TERMAN
         READ (7) RSFC
         READ (7) SEC
         READ (7) AMINUT
         READ (7) RDT
         READ (7) N
         READ (7) ISEC
         READ (7) NRAN
         READ (7) RI180
         READ (7) RIOLD180
         READ (7) SUW
         READ (7) SVW
         READ (7) SWT
         READ (7) SWQ

      ELSEIF (FLAG .EQ. 'WRITE') THEN

         WRITE (3) DPT
         WRITE (3) DQV
         WRITE (3) QCL
         WRITE (3) QRN
         WRITE (3) QCI
         WRITE (3) QCS
         WRITE (3) QCG
         WRITE (3) U
         WRITE (3) V
         WRITE (3) W
         WRITE (3) AK
         WRITE (3) DPT1
         WRITE (3) DQV1
         WRITE (3) QCL1
         WRITE (3) QRN1
         WRITE (3) QCI1
         WRITE (3) QCS1
         WRITE (3) QCG1
         WRITE (3) UU1
         WRITE (3) VV1
         WRITE (3) WW1
         WRITE (3) AK1
         WRITE (3) UMD
         WRITE (3) VMD
         WRITE (3) WMD
         WRITE (3) RSW
         WRITE (3) RLW
         WRITE (3) UB
         WRITE (3) VB
         WRITE (3) UB1
         WRITE (3) VB1
         WRITE (3) TB
         WRITE (3) QB
         WRITE (3) RHO
         WRITE (3) RHO1
         WRITE (3) TA
         WRITE (3) QA
         WRITE (3) TA1
         WRITE (3) QA1
         WRITE (3) PI
         WRITE (3) P0
         WRITE (3) WB
         WRITE (3) WBT
         WRITE (3) UBT
         WRITE (3) VBT
         WRITE (3) FD
         WRITE (3) FE
         WRITE (3) SQC
         WRITE (3) SQR
         WRITE (3) SQI
         WRITE (3) SQS
         WRITE (3) SQG
         WRITE (3) STQC
         WRITE (3) STQR
         WRITE (3) STQI
         WRITE (3) STQS
         WRITE (3) STQG
         WRITE (3) Q1T
         WRITE (3) Q2T
         WRITE (3) AQ1T
         WRITE (3) AQ2T
         WRITE (3) AQ1ZT
         WRITE (3) AQ2ZT
         WRITE (3) TAIRSFC
         WRITE (3) QAIRSFC
         WRITE (3) THAIRSF
         WRITE (3) PAIRSFC
         WRITE (3) TSDT
         WRITE (3) QSDT
         WRITE (3) THSDT
         WRITE (3) PSDT
         WRITE (3) ISUN
         WRITE (3) RLAT
         WRITE (3) RMONTH
         WRITE (3) RIDAY
         WRITE (3) IRICONT
         WRITE (3) HRL
         WRITE (3) SO0
         WRITE (3) COSZ
         WRITE (3) ICOSZ
         WRITE (3) TERMAN
         WRITE (3) RSFC
         WRITE (3) SEC
         WRITE (3) AMINUT
         WRITE (3) RDT
         WRITE (3) N
         WRITE (3) ISEC
         WRITE (3) NRAN
         WRITE (3) RI180
         WRITE (3) RIOLD180
         WRITE (3) SUW
         WRITE (3) SVW
         WRITE (3) SWT
         WRITE (3) SWQ

      ENDIF

      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE WTAPS (X,N,ifile)
      DIMENSION X(N)
      save
CC    **********************
      WRITE(ifile) X
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WCLOUD (N,AMINUT,icot,jcot,CCON,SMT,SMV,SMQ,SMC,SME,
     1   SMR,SMQC,SMQR,SMQI,SMQS,SMQG,sddd,ssss,shhh,sccc,smmm,sfff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (NX=130,NY=130,NZ=43,NT=2880,NT2=2*NT)
      COMMON/SYMB/ BLK,DOT,COM,EXM,PICE,AST,AGO,ABT
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFT,LIN,IRF,
     1   IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/BSTS5/ ACO5,ACO15,AAN5,AAN15,ANV(NT2),CNV(NT2),SPN(NT2)
      COMMON/BSTS6/ LCONV5,LANVL5,LNSPT5,IVV(NX,NY),ICS5(NX,NY,4),
     1  IBZ(NX,NY,4)
      COMMON/BI/ IT(NX,NY),ICS(NX,NY,4)
      save
      character*1 xy(nx,ny)
      character*1 blk,dot,com,exm,pice,ast,ago,abt  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(6,101) N,AMINUT,SMT,SMV,SMQ,SMC,SME,SMR
      WRITE(6,102) N,AMINUT,SMQC,SMQR,SMQI,SMQS,SMQG
      write(6,107) n,aminut,sddd,ssss,shhh,sccc,smmm,sfff
      WRITE(6,106) ACO5,AAN5,ACO15,AAN15,LCONV5,LANVL5,LNSPT5
c       ICOT=2
c       jcot=2
c$doacross local(j,i)
      DO J=2,JLES
      DO I=2,ILES
         XY(I,J)=blk
      ENDDO
      ENDDO
      DO 10 J1=2,JLES,JCOT
        J=JLES+2-J1
       WRITE(6,105) (IVV(I,J),I=2,iles,ICOT)
   10 CONTINUE
       KWATERT=11
       IF(IWATER .EQ. 1) KWATERT=18
       KCOTW=3
       IF(IWATER .EQ. 1) KCOTW=4
      DO 200 K=2,KWATERT,KCOTW
       KM=K-1
       WRITE(6,203) KM
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
       XY(I,J)=BLK
        IF(QCL(I,J,K).GE.CCON) XY(I,J)=DOT
        IF(QRN(I,J,K).GE.CCON) XY(I,J)=EXM
        IF(QRN(I,J,K).GE.CCON.AND.QCL(I,J,K).GE.CCON) XY(I,J)=COM
       ENDDO
       ENDDO
      DO J1=2,JLES,jcot
        J=JLES+2-J1
        WRITE(6,204) (XY(I,J),I=2,iles,ICOT)
      ENDDO
  200 CONTINUE
      if (iwater .eq. 0) then
        DO 400 K=12,18,6
          KM=K-1
         WRITE(6,203) KM
c$doacross local(j,i)
         DO J=2,JLES
         DO I=2,ILES
          XY(I,J)=BLK
           IF(QCI(I,J,K).GE.CCON) XY(I,J)=PICE
           IF(QCG(I,J,K).GE.CCON) XY(I,J)=AGO
           IF(QCG(I,J,K).GE.CCON.AND.QCI(I,J,K).GE.CCON) XY(I,J)=ABT
         ENDDO
         ENDDO
         DO J1=2,JLES,jcot
           J=JLES+2-J1
          WRITE(6,204) (XY(I,J),I=2,iles,ICOT)
         ENDDO
  400   CONTINUE
        DO 600 K=10,16,6
          KM=K-1
         WRITE(6,203) KM
c$doacross local(j,i)
         DO J=2,JLES
         DO I=2,ILES
          XY(I,J)=BLK
          IF(QCS(I,J,K).GE.CCON) XY(I,J)=AST
          IF(QCG(I,J,K).GE.CCON) XY(I,J)=AGO
          IF(QCS(I,J,K).GE.CCON.AND.QCG(I,J,K).GE.CCON) XY(I,J)=ABT
         ENDDO
         ENDDO
         DO J1=2,JLES,jcot
          J=JLES+2-J1
          WRITE(6,204) (XY(I,J),I=2,iles,ICOT)
         ENDDO
  600   CONTINUE
      endif
  101 FORMAT(' TIME STEP',I5,' TIME=',F8.2,' MIN SMT=',E11.4,' SMV=',
     1  E11.4,' SMQ=',E11.4,' SMC=',E11.4,' SME=',E11.4,' SMR=',E11.4)
  102 FORMAT(' TIME STEP',I5,' TIME=',F8.2,' MIN SQC=',E11.4,' SQR=',
     1  E11.4,' SQI=',E11.4,' SQS=',E11.4,' SQG',E11.4)
  105 FORMAT(1X,128I1)
  106 FORMAT(1X,' CON=',1PE10.4,1X,' ANV=',1PE10.4,1X,' CON=',1PE10.4,
     1 1X,' ANV=',1PE10.4,' CON=',I5,1X,' ANV=',I5,1X,' NONSF=',I5)
  107 format(' time step',i5,' time=',f8.2,' min dep=',e11.4,' sub=',
     1  e11.4,' qrs=',e11.4,' qrl=',e11.4,' mlt=',e11.4,' frz=',e11.4)
  203 FORMAT(24X,'K LEVEL =',I6)
  204 FORMAT(1X,128A1)
      RETURN
      END
      SUBROUTINE WCLOUDW (N,AMINUT,icot,jcot)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880)
      PARAMETER (NM16=16*NM,NT2=2*NT)
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFT,LIN,IRF,
     1   IADVH,IRFG,ISMG,ID
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5)
      COMMON/B1CR/ QCL(NX,NY,NZ),QRN(NX,NY,NZ)
      COMMON/B1IG/ QCI(NX,NY,NZ),QCG(NX,NY,NZ)
      COMMON/B1S/ QCS(NX,NY,NZ)
      COMMON/B4WP/ WW1(NX,NY,NZ)
      COMMON/BQ/ RI(NX,NY),AR(NX,NY),RX(NX,NY),RWMAX(NT2),RWMIN(NT2)
      COMMON/BBA/ XY(NX,NY),XYP(NX,NY),XYM(NX,NY),XYO(NX,NY),Y1(NM16)
      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       DO J=2,JLES
       DO I=2,ILES
        BIG=0.
        SMALL=0.
        QMAXL=0.
        QMAXI=0.
        DO 10 K=2,KLES
          XA=WW1(I,J,K)
         BIG=MAX(BIG,XA)
         SMALL=MIN(SMALL,XA)
          XB=QCL(I,J,K)+QRN(I,J,K)
         QMAXL=MAX(QMAXL, XB)
          XC=QCI(I,J,K)+QCS(I,J,K)+QCG(I,J,K)
         QMAXI=MAX(XC, QMAXI)
   10   CONTINUE
       XY(I,J)=.01*BIG
       XYP(I,J)=.01*ABS(SMALL)
       XYM(I,J)=1.E3*QMAXL
       XYO(I,J)=1.E3*QMAXI
      ENDDO
      ENDDO

       WRITE(6,101) N,AMINUT
       CALL WOUT (RI,icot,jcot,1)
       WRITE(6,102) N,AMINUT
       CALL WOUT (XY,icot,jcot,2)
       WRITE(6,103) N,AMINUT
       CALL WOUT (XYP,icot,jcot,3)
       WRITE(6,104) N,AMINUT
       CALL WOUT (XYM,icot,jcot,4)
       if (iwater .eq. 0) then
         WRITE(6,105) N,AMINUT
         CALL WOUT (XYO,icot,jcot,5)
       endif
  101 FORMAT(' TIME STEP',I5,' TIME=',F8.2,14X,'RAINFALL')
  102 FORMAT(' TIME STEP',I5,' TIME=',F8.2,14X,'UPDRAFT')
  103 FORMAT(' TIME STEP',I5,' TIME=',F8.2,14X,'DOWNDRAFT')
  104 FORMAT(' TIME STEP',I5,' TIME=',F8.2,14X,'LIQUID WATER')
  105 FORMAT(' TIME STEP',I5,' TIME=',F8.2,14X,'ICE WATER')
      RETURN
      END

      SUBROUTINE WOUT (XY,icot,jcot,IDD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (NX=130,NY=130)
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,NISEC(5) 
      DIMENSION XY(NX,NY),iXY(NX,NY),IWH(10)
      save
      DATA IWH/1,2,3,4,5,6,7,8,9,0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       A1=0.01
       A2=2.
       A3=4.
       A4=6.
       A5=8.
       A6=10.
       A7=12.
       A8=14.
       A9=16.
       A10=18.
       IF(IDD.EQ.1) THEN
        A2=1.
        A3=5.
        A4=10.
        A5=15.
        A6=20.
        A7=25.
        A8=30.
        A9=35.
        A10=40.
       ENDIF
       IF(IDD.EQ.3) THEN
        A2=0.5
        A3=1.0
        A4=1.5
        A5=2.0
        A6=2.5
        A7=3.
        A8=3.5
        A9=4.
        A10=5.
       ENDIF
       IF(IDD.GE.4) THEN
        A2=0.1
        A3=0.5
        A4=1.
        A5=2.
        A6=3.
        A7=4.
        A8=5.
        A9=7.5
        A10=10.
       ENDIF
c$doacross local(j,i)
       DO J=2,JLES
       DO I=2,ILES
       IXY(I,J)=IWH(10)
        IF(XY(I,J).GE.A1) IXY(I,J)=IWH(1)
        IF(XY(I,J).GE.A2) IXY(I,J)=IWH(2)
        IF(XY(I,J).GE.A3) IXY(I,J)=IWH(3)
        IF(XY(I,J).GE.A4) IXY(I,J)=IWH(4)
        IF(XY(I,J).GE.A5) IXY(I,J)=IWH(5)
        IF(XY(I,J).GE.A6) IXY(I,J)=IWH(6)
        IF(XY(I,J).GE.A7) IXY(I,J)=IWH(7)
        IF(XY(I,J).GE.A8) IXY(I,J)=IWH(8)
        IF(XY(I,J).GE.A9) IXY(I,J)=IWH(9)
        IF(XY(I,J).GE.A10) IXY(I,J)=IWH(9)
       ENDDO
       ENDDO
      DO J1=2,JLES,JCOT
        J=JLES+2-J1
        WRITE(6,204) (IXY(I,J),I=2,iles,ICOT)
      ENDDO
  204 FORMAT(1X,128I1)
      RETURN
      END

C
      SUBROUTINE RINIT (IRS)
C     ******   INITIALIZE ALL VARIABLES TO 0.0S   ******
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880,ITT=244)
      PARAMETER (NT2=2*NT,NB0=2*NZ+14,NB1=178,NB2=193)
      PARAMETER (NB3=5*NZ+8*NX*NY,NB4=21*NZ+NX,NB5=11*NZ,NB6=5*NZ)
      PARAMETER (NB7=13*NZ,NB8=62*NZ,NB10=12*NZ+5*NT)
      PARAMETER (NB11=4*NX*NY+16*NM)
      PARAMETER (NB14=3*NX*NY+2*NT2,NB16=4*NZ+1)
      PARAMETER (NB17=48*NZ*4*7+2*NZ*4,NB18=13*NZ*4*7,NB19=12*NZ*4*7)
      PARAMETER (NB20=4*NZ*4*7,NB21=4+3*NT2,NB22=3+NX*NY+2*NX*NY*4)
      PARAMETER (NB23=NX*NY+NX*NY*4,NB24=NZ*16*14+NZ*5*14)
      PARAMETER (NB25=16*NZ*14+5*NZ*14+2*NZ,NB26=16*NZ*14+5*NZ*14)
      parameter (nb27=2*nz*4*7,nxy2=2*nx*ny,nz6=6*nz)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      parameter (nz12=12*nz,nb66=5*itt*nz+5*nz,nb666=2*nz*itt+2*nz,
     1           nz9=9*nz)
      common/bb6/ tls(nz12)  
      common/bb66/ u66h(nb66)
      COMMON/q1q2z/ Q666H(Nb666)
      common/gbs/ tlsw(nz9)
      COMMON/BXYZ/ IMAX(14)
      COMMON/BX/ VX(9)
      COMMON/BY/ VY(9)
      COMMON/BZ/ VZ(9)
      COMMON/CONT/ RCONT(11)
      COMMON/BT/ RDT(19)
      COMMON/O4X/ OX4(NB0)
      COMMON/BTERV/ VTERV(6)
      COMMON/BSNW/ VSNW(NB1)
      COMMON/RTERV/ RTER(9)
      COMMON/RSNW/ RSN(NB2)
CC
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
cc      COMMON/BW/ TRAHX(NX,NY,NZ),TRAHY(NX,NY,NZ),TRAV(NX,NY,NZ)
      COMMON/BSAT/ XYZ1(NX,NY,NZ)
      COMMON/BSAT1/ XYZ2(NX,NY,NZ)
      COMMON/BADV/ XYZ3(NX,NY,NZ)
      COMMON/BADV1/ XYZ4(NX,NY,NZ)
      COMMON/B1TQ/ V11(NX,NY,NZ),V12(NX,NY,NZ)
      COMMON/B1CR/ V13(NX,NY,NZ),V14(NX,NY,NZ)
      COMMON/B1IG/ V15(NX,NY,NZ),V16(NX,NY,NZ)
      COMMON/B2TQ/ VB1(NX,NY,NZ),VB2(NX,NY,NZ)
      COMMON/B2CR/ VB3(NX,NY,NZ),VB4(NX,NY,NZ)
      COMMON/B2IG/ VB5(NX,NY,NZ),VB6(NX,NY,NZ)
      COMMON/B3UV/ V21(NX,NY,NZ),V22(NX,NY,NZ)
      COMMON/B4UV/ VB21(NX,NY,NZ),VB22(NX,NY,NZ)
      COMMON/B1S/ V17(NX,NY,NZ)
      COMMON/B2S/ VB7(NX,NY,NZ)
      COMMON/B3WP/ V23(NX,NY,NZ)
      COMMON/B4WP/ VB23(NX,NY,NZ)
      COMMON/B3A/ V41(NX,NY,NZ)
      COMMON/B4A/ V42(NX,NY,NZ)
      COMMON/B4/ TBSK(NB3)
      COMMON/B5/ TB(NB4)
      COMMON/B6/ FD(NB5)
      COMMON/B7/ SQAAQ(NB6)
      COMMON/B8/ SMF(NB7)
      COMMON/B9/ SUT1(NB8)
      COMMON/BRH1/ SRRO(NB10)
      COMMON/BBA/ XYP(NB11)
      COMMON/BQ/ RI(NB14)
      COMMON/DAMP/ RFA(NB16)
      COMMON/BSKD/ SKCOEF(NZ)

      COMMON/BSTS/ STS(NB17)
      COMMON/BSTS1/ STS1(NB18)
      COMMON/BSTS2/ STS2(NB18)
      COMMON/BSTS3/ STS3(NB19)
      COMMON/BSTS4/ STS4(NB20)
      COMMON/BSTS5/ STS5(NB21)
      COMMON/BSTS6/ LBST(NB22)
      COMMON/BSTS20/ STS20(NB27)
      COMMON/BI/ IT(NB23)
      COMMON/BTT/ S1(NB24)
      COMMON/BCS/ S2(NB25)
      COMMON/BCSS/ S22(NB26)

      COMMON/BSFC/ tsfc_1(nxy2)
      common/b66b/ s_dep(nz6)
      common/TMP/ y1d(nx,ny,nz),y2d(nx,ny,nz),y3d(nx,ny,nz),
     1            y4d(nx,ny,nz),y03d(nx,ny,1)
      COMMON /tmp1/ U1(NX,NY,NZ),V1(NX,NY,NZ),W1(NX,NY,NZ)

C_TAO

      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$doacross local(k,j,i)
      DO K=1,NZ
         DO J=1,NY
            DO I=1,NX
               Y1D(i,j,k)=0.
               Y2D(i,j,k)=0.
               Y3D(i,j,k)=0.
               Y4D(i,j,k)=0.
               U1(i,j,k)=0.
               V1(i,j,k)=0.
               W1(i,j,k)=0.
            ENDDO
         ENDDO
      ENDDO

c$doacross local(j,i)
      DO J=1,NY
        DO I=1,NX
            y03d(i,j,1)=0.
         ENDDO
      ENDDO

      do i=1,nb27
        sts20(i)=0.
      enddo
      DO I=1,NB17
         STS(I)=0.0
      ENDDO
      DO I=1,NB18
        STS1(I)=0.0
        STS2(I)=0.0
      ENDDO
      DO I=1,NB19
        STS3(I)=0.0
      ENDDO
      DO I=1,NB20
        STS4(I)=0.0
      ENDDO
      DO I=1,NB21
        STS5(I)=0.0
      ENDDO
      DO I=1,NB22
        LBST(I)=0
      ENDDO
      DO I=1,NB23
         IT(I)=0
      ENDDO
      DO I=1,NB24
        S1(I)=0.
      ENDDO
      DO I=1,NB25
        S2(I)=0.
      ENDDO
      DO I=1,NB26
         S22(I)=0.
      ENDDO
      do k=1,nz12
         tls(k)=0.
      enddo
      do k=1,nb66
         u66h(k)=0.
      enddo
      do k=1,nb666
         Q666H(k)=0.
      enddo
      do k=1,nz9
         tlsw(k)=0.
      enddo
      do i=1,nxy2
         tsfc_1(i)=0.
      enddo
       IF(IRS.EQ.1) RETURN
      DO I=1,14
        IMAX(I)=0
      ENDDO
      DO I=1,9
       RTER(I)=0.
       VX(I)=0.
       VY(I)=0.
       VZ(I)=0.
      ENDDO
      DO I=1,11
        RCONT(I)=0.0
      ENDDO
      DO I=1,19
         RDT(I)=0.0
      ENDDO
      DO I=1,NB0
         OX4(I)=0.0
      ENDDO
      DO I=1,6
        VTERV(I)=0.
      ENDDO
      DO I=1,NB1
        VSNW(I)=0.
      ENDDO
      DO I=1,NB2
         RSN(I)=0.
      ENDDO
      DO 1002 K=1,NZ
      DO 1001 J=1,NY
      DO 1000 I=1,NX

       UMD(I,J,K)=0.
       VMD(I,J,K)=0.
       WMD(I,J,K)=0.
c       TRAHX(I,J,K)=0.
c       TRAHY(I,J,K)=0.
c       TRAV(I,J,K)=0.
       XYZ1(I,J,K)=0.0
       XYZ2(I,J,K)=0.0
       XYZ3(I,J,K)=0.0
       XYZ4(I,J,K)=0.0
       V11(I,J,K)=0.0
       V12(I,J,K)=0.0
       V13(I,J,K)=0.0
       V14(I,J,K)=0.0
       V15(I,J,K)=0.0
       V16(I,J,K)=0.0
       V17(I,J,K)=0.0
       VB1(I,J,K)=0.0
       VB2(I,J,K)=0.0
       VB3(I,J,K)=0.0
       VB4(I,J,K)=0.0
       VB5(I,J,K)=0.0
       VB6(I,J,K)=0.0
       VB7(I,J,K)=0.0
       V21(I,J,K)=0.0
       V22(I,J,K)=0.0
       V23(I,J,K)=0.0
       VB21(I,J,K)=0.0
       VB22(I,J,K)=0.0
       VB23(I,J,K)=0.0
       V41(I,J,K)=0.0
       V42(I,J,K)=0.0
 1000 CONTINUE
 1001 CONTINUE
 1002 CONTINUE
      DO I=1,NB3
        TBSK(I)=0.0
      ENDDO
      DO I=1,NB4
        TB(I)=0.0
      ENDDO
      DO I=1,NB5
        FD(I)=0.0
      ENDDO
      DO  I=1,NB6
        SQAAQ(I)=0.0
      ENDDO
      DO I=1,NB7
        SMF(I)=0.0
      ENDDO
      DO I=1,NB8
        SUT1(I)=0.0
      ENDDO
      DO I=1,NB10
        SRRO(I)=0.0
      ENDDO
      DO I=1,NB11
        XYP(I)=0.
      ENDDO
      DO I=1,NB14
        RI(I)=0.
      ENDDO
      DO I=1,NB16
        RFA(I)=0.
      ENDDO
      DO K=1,NZ
        SKCOEF(K)=0.
      ENDDO
      do k=1,nz6
         s_dep(k)=0.
      enddo
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TERP1 (N,X,F,W,Y,INT,TAB,ITAB)
      DIMENSION X(n),F(n),W(n),TAB(3),ITAB(3)
      save
      CALL SEARCH (Y,X,N,I)
      CALL INTERP (N,X,F,W,Y,I,INT,TAB,ITAB)
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SEARCH (XBAR,X,N,I)
      DIMENSION X(n)
      save
      DATA B/.69314718/
      IF(XBAR.GT.X(2)) GO TO 101
      I=1
      RETURN
  101 CONTINUE
      IF(XBAR.LT.X(N-1)) GO TO 102
      I=N-1
      RETURN
  102 CONTINUE
      M=INT((LOG(FLOAT(N)))/B)
      I=2**M
      IF(I.GE.N) I=I/2
      K=I
      NM1=N-1
  103 CONTINUE
      K=K/2
      IF(XBAR.GE.X(I)) GO TO 104
      I=I-K
      GO TO 103
  104 CONTINUE
      IF(XBAR.LE.X(I+1)) RETURN
      I=MIN0(I+K,NM1)
      GO TO 103
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE INTERP (N,X,F,W,Y,I,INT,TAB,ITAB)
      DIMENSION X(n),F(n),W(n),TAB(3),ITAB(3)
      save
      II(I)=(I-1)*INT+1
      FLK=X(I+1)-X(I)
      FLP=X(I+1)-Y
      FL0=Y-X(I)
      I0=II(I)
      IP=I0+INT
      IF(ITAB(1).NE.1) GO TO 102
      A=(W(I0)*FLP**3+W(IP)*FL0**3)/(6.*FLK)
      B=(F(IP)/FLK-W(IP)*FLK/6.)*FL0
      C=(F(I0)/FLK-W(I0)*FLK/6.)*FLP
      TAB(1)=A+B+C
  102 IF(ITAB(2).NE.1) GO TO 104
      A=(W(IP)*FL0**2-W(I0)*FLP**2)/(2.*FLK)
      B=(F(IP)-F(I0))/FLK
      C=(W(I0)-W(IP))*FLK/6.
      TAB(2)=A+B+C
  104 IF(ITAB(3).NE.1) GO TO 106
      TAB(3)=(W(I0)*FLP+W(IP)*FL0)/FLK
  106 RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE COEFF (N,X,F,W,IOP,INT,WK)
      DIMENSION X(n),F(n),W(n),IOP(2),WK(N,4)
      save
      II(I)=(I-1)*INT+1
      J0=1
      DO 101 I=2,N
      JM=J0
      J0=J0+INT
      WK(I,1)=X(I)-X(I-1)
      WK(I,2)=(F(J0)-F(JM))/WK(I,1)
      WK(I,3)=WK(I,1)/6.
      WK(I,1)=WK(I,1)/3.
  101 CONTINUE
      NN=N
      MK=IOP(1)
      ML=IOP(2)
      GO TO (102,103,104,105) ,MK
  102 CONTINUE
      WK(2,2)=WK(3,2)-WK(2,2)-WK(2,3)*W(1)
      WK(2,3)=0.
      WK(2,1)=WK(2,1)+WK(3,1)
      I1=2
      NN=NN-1
      GO TO 106
  103 CONTINUE
      WK(1,2)=WK(2,2)-W(1)
      WK(2,2)=WK(3,2)-WK(2,2)
      WK(1,3)=0.
      WK(1,1)=WK(2,1)
      WK(2,1)=WK(2,1)+WK(3,1)
      I1=1
      GO TO 106
  104 CONTINUE
      Y2=WK(2,2)
      B2=WK(2,1)
      WK(2,2)=WK(3,2)-WK(2,2)
      WK(2,1)=WK(3,1)+WK(2,1)
      I1=2
      NN=NN-1
      GO TO 106
  105 CONTINUE
      A12=X(1)-X(2)
      A13=X(1)-X(3)
      A14=X(1)-X(4)
      A23=X(2)-X(3)
      A24=X(2)-X(4)
      A34=X(3)-X(4)
      J1=1
      J2=J1+INT
      J3=J2+INT
      J4=J3+INT
      W(1)=(1./A12+1./A13+1./A14)*F(J1)-
     1     A13*A14/(A12*A23*A24)*F(J2)+A12*A14/(A13*A23*A34)*F(J3)-
     2     A12*A13/(A14*A24*A34)*F(J4)
      GO TO 103
  106 CONTINUE
      I2=N-2
      DO 107 I=3,I2
      WK(I,2)=WK(I+1,2)-WK(I,2)
      WK(I,1)=WK(I+1,1)+WK(I,1)
  107 CONTINUE
      IN=II(N)
      GO TO (108,109,110,111) ,ML
  108 CONTINUE
      WK(N-1,2)=WK(N,2)-WK(N-1,2)-WK(N,3)*W(IN)
      WK(N,3)=0.
      WK(N-1,1)=WK(N-1,1)+WK(N,1)
      NN=NN-1
      GO TO 112
  109 CONTINUE
      WK(N-1,2)=WK(N,2)-WK(N-1,2)
      WK(N,2)=-WK(N,2)+W(IN)
      WK(N-1,1)=WK(N-1,1)+WK(N,1)
      WK(1,4)=0.
      GO TO 112
  110 CONTINUE
      WK(N-1,2)=WK(N,2)-WK(N-1,2)
      WK(N,2)=Y2-WK(N,2)
      WK(N-1,1)=WK(N-1,1)+WK(N,1)
      WK(N,1)=WK(N,1)+B2
      WK(1,4)=WK(2,3)
      GO TO 112
  111 CONTINUE
      A12=X(N)-X(N-1)
      A13=X(N)-X(N-2)
      A14=X(N)-X(N-3)
      A23=X(N-1)-X(N-2)
      A24=X(N-1)-X(N-3)
      A34=X(N-2)-X(N-3)
      J1=IN
      J2=J1-INT
      J3=J2-INT
      J4=J3-INT
      W(IN)=(1./A12+1./A13+1./A14)*F(J1)-
     1      A13*A14/(A12*A23*A24)*F(J2)+A12*A14/(A13*A23*A34)*F(J3)-
     2      A12*A13/(A14*A24*A34)*F(J4)
      GO TO 109
  112 CONTINUE
      II1=II(I1)
      CALL TRIP (NN,WK(I1,3),WK(I1,1),WK(I1+1,3),WK(I1,2),W(II1),INT)
      GO TO (114,114,113,114) ,MK
  113 CONTINUE
      W(1)=W(IN)
  114 CONTINUE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TRIP (N,A,B,C,Y,Z,INT)
      DIMENSION A(n),B(n),C(n),Y(n),Z(n)
      save
      II(I)=(I-1)*INT+1
      BN=B(N)
      YN=Y(N)
      V=C(N)
      Y(1)=Y(1)/B(1)
      A(1)=A(1)/B(1)
      B(1)=C(1)/B(1)
      NM2=N-2
      DO 101 J=2,NM2
      DEN=B(J)-A(J)*B(J-1)
      B(J)=C(J)/DEN
      Y(J)=(Y(J)-A(J)*Y(J-1))/DEN
      A(J)=-A(J)*A(J-1)/DEN
      BN=BN-V*A(J-1)
      YN=YN-V*Y(J-1)
      V=-V*B(J-1)
  101 CONTINUE
      DEN=B(N-1)-A(N-1)*B(N-2)
      B(N-1)=(C(N-1)-A(N-1)*A(N-2))/DEN
      Y(N-1)=(Y(N-1)-A(N-1)*Y(N-2))/DEN
      BN=BN-V*A(N-2)
      YN=YN-V*Y(N-2)
      V=A(N)-V*B(N-2)
      NM1=N-1
      IN=II(N)
      INM=II(NM1)
      Z(IN)=(YN-V*Y(NM1))/(BN-V*B(NM1))
      Z(INM)=Y(NM1)-B(NM1)*Z(IN)
      DO 102 J=2,NM1
      K=N-J
      IK=II(K)
      Z(IK)=Y(K)-B(K)*Z(IK+INT)-A(K)*Z(IN)
  102 CONTINUE
      RETURN
      END


      SUBROUTINE CONSATRH
C    (R&H)  SPECIFY SOME CONSTANTS IN SATICE ROUTINE   ******
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK(8)
      COMMON/SIZE/ TNW,TNS,TNG,ROQS,ROQG,ROQR
      COMMON/RTERV/ ZRC,ZGC,ZSC,VRC0,VRC1,VRC2,VRC3,VGC,VSC
      COMMON/B3CS/ AG,BG,AS,BS,AW,BW,BGH,BGQ,BSH,BSQ,BWH,BWQ
      COMMON/RSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,RN2,BND2,RN3,RN4,
     1  RN5,RN50,RN51,RN52,RN53,RN6,RN60,RN61,RN62,RN63,RN7,RN8,RN9,
     2  RN10,RN101,RN102,RN10A,RN10B,RN10C,RN11,RN12,RN12A(31),
     3  RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN171,RN172,RN17A,RN17B,
     4  RN17C,RN18,RN18A,RN19,RN191,RN192,RN19A,RN20,RN20A,RN20B,RN30,
     5  RN30A,RN21,BND21,RN22,RN23,RN231,RN232,RN25,RN25A(31),RN31,BETA,
     6  RN32,RN33,RN331,RN332,RN34,RN35
      DIMENSION A1(31),A2(31)
      save
      DATA A1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,
     1   .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,
     2   .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,
     3   .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,
     4   .6513E-6,.5956E-6,.5333E-6,.4834E-6/
      DATA A2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
C     ******************************************************************
        CPI=4.*ATAN(1.)
        CPI2=CPI*CPI
c        GRVT=980.
      C38=3.799052E3
      C358=35.86
      C610=6.1078E3
      C149=1.496286E-5
      C879=8.794142
      C172=17.26939
      C409=4098.026
      C76=7.66
      C218=21.87456
      C580=5807.695
      C141=1.414435E7
C     ***************
         TCA=2.43E3
         DWV=.226
         DVA=1.718E-4
         AMW=18.016
         ARS=8.314E7
      T0=273.16
      T00=238.16
      ALV=2.5E10
      ALF=3.336E9
      ALS=2.8336E10
       AVC=ALV/CP
       AFC=ALF/CP
       ASC=ALS/CP
      RW=4.615E6
      CW=4.187E7
      CI=2.093E7
C***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
      ROQR=1.
      TNW=.08
        ROQS=.1
        TNS=.16
CFRED   TNS=1.
          ROQG=.4
           TNG=.08
C        TNG=.22
C***   DEFINE THE COEFFICIENTS USED IN TERMINAL VELOCITY
       AG=351.2
       BG=.37
c        AS=132.93
c        BS=.11
c         AS=152.93
c        AS=159.
c        BS=.25
          as=78.63154
          bs=.11

         AW=2115.
         BW=.8
       BGH=.5*BG
       BSH=.5*BS
       BWH=.5*BW
       BGQ=.25*BG
       BSQ=.25*BS
       BWQ=.25*BW
      GA3=2.
      GA4=6.
      GA5=24.
      GA6=120.
      GA7=720.
      GA8=5040.
      GA9=40320.
        GA4G=11.63177
        GA3G=3.3233625
        GA5GH=1.608355
        IF(BG.EQ.0.37) GA4G=9.730877
        IF(BG.EQ.0.37) GA3G=2.8875
        IF(BG.EQ.0.37) GA5GH=1.526425
          GA3D=2.54925
          GA4D=8.285063
          GA5DH=1.456943
          IF(BS.EQ.0.57) GA3D=3.59304
          IF(BS.EQ.0.57) GA4D=12.82715
          IF(BS.EQ.0.57) GA5DH=1.655588
          IF(BS.EQ.0.11) GA3D=2.218906
          IF(BS.EQ.0.11) GA4D=6.900796
          IF(BS.EQ.0.11) GA5DH=1.382792
          GA6D=144.93124
CCCCCC             RUTLEDGE AND HOBBS, 1984   CCCCCCCCCCCCCCCCCCCCCCCCCC
       AC1=AS
       AC2=AG
       ZRC=(CPI*ROQR*TNW)**0.25
       ZSC=(CPI*ROQS*TNS)**0.25
       ZGC=(CPI*ROQG*TNG)**0.25
       VRC0=-26.7
       VRC1=20600./ZRC
       VRC2=-204500./(ZRC*ZRC)
       VRC3=906000./(ZRC*ZRC*ZRC)
ctao       VSC=1.25*AC1*GA4D/(6.*ZSC**BS)
       VSC=AC1*GA4D/(6.*ZSC**BS)
       VGC=AC2*GA4G/(6.*ZGC**BG)
CS      CD1=6.E-1
CS      CD2=4.*GRVT/(3.*CD1)
CS     VGC=GA4G*SQRT(CD2*ROQG/ZGC)/6.
C     ****************************
      RN1=9.4E-15
      RN2=1.E-3
       BND2=2.E-3
       ESI=.1
      RN3=.25*CPI*TNS*AC1*ESI*GA3D
       ESC=1.
      RN4=.25*CPI*ESC*TNS*AC1*GA3D
CRH    ERI=1.
       ERI=.1   
      RN5=.25*CPI*ERI*TNW
       RN50=-.267E2*GA3
       RN51=5.15E3*GA4
       RN52=-1.0225E4*GA5
       RN53=7.55E3*GA6
       AMI=1./(24.*6.E-9)
      RN6=CPI2*ERI*TNW*ROQR*AMI
       RN60=-.267E2*GA6
       RN61=5.15E3*GA7
       RN62=-1.0225E4*GA8
       RN63=7.55E3*GA9
       ESR=1.
       esr=0.5*esr
      RN7=CPI2*ESR*TNW*TNS*ROQS
       esr=1.
      RN8=CPI2*ESR*TNW*TNS*ROQR
       EGS=.1
      RN9=CPI2*EGS*TNS*TNG*ROQS
      RN10=4.*TNS
       RN101=.65
       RN102=.44*SQRT(AC1/DVA)*GA5DH
       RN10A=ALV*ALS*AMW/(TCA*ARS)
       RN10B=ALV/TCA
       RN10C=ARS/(DWV*AMW)
      RN11=2.*CPI*TNS*TCA/ALF
CLIN   AMI50=4.8E-7
       AMI50=4.8E-7*(100./50.)**3
       AMI40=2.46E-7
       ami40=2.46e-7*.5**3
       EIW=1.
       UI50=100.
CLIN   RI50=5.E-3
       RI50=(100./50.)*5.E-3
       CMN=1.05E-15
      RN12=CPI*EIW*UI50*RI50*RI50
      DO 10 K=1,31
        Y1=1.-A2(K)
       RN13(K)=A1(K)*Y1/(AMI50**Y1-AMI40**Y1)
        RN12A(K)=RN13(K)/AMI50
        RN12B(K)=A1(K)*AMI50**A2(K)
        RN25A(K)=A1(K)*CMN**A2(K)
   10 CONTINUE
       EGC=1.
      RN14=.25*CPI*EGC*AC2*TNG*GA3G
       EGI=.1
      RN15=.25*CPI*EGI*TNG*AC2*GA3G
       EGI=1.
      RN15A=.25*CPI*EGI*TNG*AC2*GA3G
       EGR=1.
      RN16=CPI2*EGR*TNG*TNW*ROQR
      RN171=2.*CPI*TNG*ALV*DWV
       RN172=2.*CPI*TNG*TCA
       RN17A=.31*GA5GH*SQRT(AC2/DVA)
       RN17B=CW-CI
       RN17C=CW
       APRI=.66
       BPRI=1.E-4
       BPRI=0.5*BPRI
      RN18=20.*CPI2*BPRI*TNW*ROQR
       RN18A=APRI
      RN19=2.*CPI*TNG*TCA/ALF
       RN191=.78
       RN192=.31*GA5GH*SQRT(AC2/DVA)
       RN19A=CW/ALF
      RN20=2.*CPI*TNG
       RN20A=ALS*ALS*AMW/(TCA*ARS)
       RN20B=ALS/TCA
      RN30=2.*CPI*TNG
       RN30A=ALV*ALV*AMW/(TCA*ARS)
      RN21=1.E-3
cc       BND21=1.E-3
       BND21=1.5E-3
       ERC=1.
      RN22=.25*CPI*ERC*TNW
      RN23=2.*CPI*TNW
       RN231=.78
       RN232=.31*GA3*SQRT(3.E3/DVA)
       CN0=1.E-8
      rn25=cn0/1000.
      RN25=CN0
      RN31=1.E-17
       BETA=-.6
      RN32=4.*51.545E-4
      RN33=4.*TNS
       RN331=.65
       RN332=.44*SQRT(AC1/DVA)*GA5DH
       ESC=1.
       ESC=0.5*ESC
       AMC=1./(24.*4.E-9)
      RN34=CPI2*ESC*AMC*AC1*ROQS*TNS*GA6D
      RN35=ALV*ALV/(CP*RW)
C     ****************************
      RETURN
      END

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION CVMGP(X1,X2,X3)
       IF(X3.GE.0.)  THEN
        CVMGP=X1
       ELSE
        CVMGP=x2
       ENDIF
      RETURN
      END
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION RAN3(idum)
      implicit none
      integer mbig,mseed,mz,idum
      real    ran3,fac
      PARAMETER  (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      integer i,k,iff,mj,mk,ii,inext,inextp,ma(55)
      SAVE  inext,inextp,ma,iff
      DATA  iff /0/
      IF( idum.LT.0 .OR. iff.EQ.0 ) THEN
        iff= 1
        mj= mseed - ABS( idum )
c       mj=MOD( mj,mbig )
        mj= mj-mbig*(mj/mbig)
        ma(55)= mj
        mk= 1
        DO 10 i=1,54
c          ii= MOD( 21*i,55 )
          ii= (21*i)-55*((21*i)/55)
          ma(ii)= mk
          mk= mj-mk
          IF( mk.LT.mz ) mk= mk+mbig
          mj= ma(ii)
 10     CONTINUE
        DO 20 k=1,4
          DO 15 i=1,55
c            ma(i)= ma(i) - ma(1+MOD( i+30,55 ))
            ma(i)= ma(i) - ma(1+(i+30)-55*((i+30)/55))
            IF( ma(i).LT.mz ) ma(i)= ma(i)+mbig
 15       CONTINUE
 20     CONTINUE
        inext= 0
        inextp= 31
        idum= 1
      ENDIF
      inext= inext+1
      IF( inext.EQ.56 ) inext= 1
      inextp= inextp+1
      IF( inextp.EQ.56 ) inextp= 1
      mj= ma(inext) - ma(inextp)
      IF( mj.LT.mz ) mj= mj+mbig
      ma(inext)= mj
      ran3= REAL( mj )*fac
      RETURN
      END
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ainter(sst,ssp,tb,fd,qb,fe,ub,vb,wb,y1,y2,topht)
      PARAMETER (N=41,ITT=244,NZ=43,KLES=NZ-1,KL=18,KMAX1=KLES+1)
      PARAMETER (KBA=18*KMAX1,NN1=N+1,NN2=2*N+1)
      common/msui/ isounding,isui,ifine
      dimension y1(nz),y2(nz),itab(3),iop(2)
      real fd(nz),fe(nz),tb(nz),qb(nz),ub(nz),vb(nz),wb(nz),uo_6h(n,itt)
     1    ,vo_6h(n,itt),wo_6h(n,itt),q1o_6h(n,itt),q2o_6h(n,itt)
      real P(N),T(N),Q(N),U(N),V(N),WW1(N),H1(N),H2(N)
      real TAB(3),rqo,AT(162),AQ(162),H(N),WK(N,4),W(N),QQ(N)
      real Z(nz),Z1(nz),ZZ(168)
c      real YYZ(180)
      real AA(KMAX1,KL),BB(KMAX1,KL),QQ1(3*N),QQ2(KBA),QQ3(KBA)
      common/bb6/ tls(nz),qls(nz),tls1(nz),tls2(nz),qls1(nz),qls2(nz),
     1   tls3(nz),tls4(nz),qls3(nz),qls4(nz),sft(nz),sfq(nz)
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1  wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2  q2t(nz)
      COMMON/UPPER/ T_ADJUST(N,ITT),Q_ADJUST(N,ITT),P_T(N),
     1      T_T(N),V_T(N)
      common/upper1/ t_adjust0(nz,itt),q_adjust0(nz,itt),press_t0(nz),
     1     temp_t0(nz),vap_t0(nz)
      common /zzzobs/ zzz(n,itt)
      real sst(itt),ssp(itt),dum(n)
      EQUIVALENCE (AA(1,1),QQ2(1)),(BB(1,1),QQ3(1))
      EQUIVALENCE (U(1),QQ1(1)),(V(1),QQ1(NN1)),(WW1(1),QQ1(NN2))
      character*4 IVAA(9)
      save
      DATA IVAA/'P  ',' H ',' T ',' Q ',' U ',' V ','WW1 ',
     1   'RHW ','RHI '/
      call obs(sst,ssp,p,t,q,u,v,ww1,uo_6h,vo_6h,wo_6h,
     1              q1o_6h,q2o_6h)
c      call obs_data(sst,ssp,p,t,q,u,v,ww1,uo_6h,vo_6h,wo_6h,q1o_6h,
c     1              q2o_6h)
      KMAX=KLES+KLES+1
c      DZ=DZZ*.01
      DO K=1,KLES
        Z(K)=y2(k+1)*.01
        Z1(K)=y1(k+1)*.01
      ENDDO
      Z(KMAX1)=Z(KLES)+(Z(KLES)-Z(KLES-1))
      Z1(KMAX1)=Z1(KLES)+(Z1(KLES)-Z1(KLES-1))
      write(6,2031)
      DO K=1,KMAX1
         write(6,100) K,Z(K),Z1(K)
      ENDDO
      DO K=1,KLES
        I=K*2
        I1=K*2-1
        ZZ(I1)=Z(K)
        ZZ(I)=Z1(K)
      ENDDO
      ZZ(KMAX)=Z(KMAX1)
      H(1)=0.
      TERM1=2.87E+6/980.*.5
      DO 20 K=2,N
      KM=K-1
c     H(K)=H(KM)-TERM1*((T(K)+273.16)*(1.+.61*Q(K)*1.E-3)
c    1     +(T(KM)+273.16)*(1.+.61*Q(KM)*1.E-3))*LOG(P(K)/P(KM))*.01
c      print*,'h,z',h(k),zzz(k,1)
      h(k)=zzz(k,1)  ! 4/9/01
       ADD=3.799E3/P(KM)*EXP(17.26939-4098.026/(T(KM)+273.16-35.86))
      H1(KM)=Q(KM)/ADD
       ADD1=3.799E3/P(KM)*EXP(21.87456-5807.695/(T(KM)+273.16-7.66))
      H2(KM)=Q(KM)/ADD1
   20 CONTINUE

      zpresst=100.
      do k=1,n
        if(p(k).ge.zpresst) topht=h(k)
      enddo
      print*,'top at or just below 100mb =',topht

c      topht=h(n-2)
      write(6,102) (IVAA(I),I=1,9)
      DO K=1,N
       write(6,1031)K,P(K),H(K),T(K),Q(K),U(K),V(K),WW1(K),H1(K),H2(K)
      ENDDO

      do k=1,n
        if(h(k).le.z(nz-1)) ntop=k
      enddo
      rqo=0.
      do k=1,ntop
        rqo=rqo+.5*(q(k)+q(k+1))*.001*(p(k)-p(k+1))*1000./980.
      enddo
       print*
       print*,'OBSERVED QV INTEGRATION=',rqo

      IOP(1)=4
      IOP(2)=4
      INT=1
      ITAB(1)=1
      ITAB(2)=0
      ITAB(3)=0
      CALL COEFF (N,H,T,W,IOP,INT,WK)
      DO K=1,KMAX
        Y=ZZ(K)
        CALL TERP1 (N,H,T,W,Y,INT,TAB,ITAB)
        AT(K)=TAB(1)
      ENDDO
      CALL COEFF (N,H,Q,W,IOP,INT,WK)
      write(6,1021) (IVAA(I),I=2,4)
      DO 250 K=1,KMAX
        Y=ZZ(K)
       CALL TERP1 (N,H,Q,W,Y,INT,TAB,ITAB)
       AQ(K)=TAB(1)
        if(y .ge. 16000.) aq(k)=0.
       write(6,106) K,Y,AT(K),AQ(K)
c        YYZ(K)=Y
  250 CONTINUE
      DO KK=1,KLES
         kkp=kk+1
        KK1=(KK-1)*2+1
        KK2=KK1+1
       FD(KKp)=AT(KK1)
       TB(KKp)=AT(KK2)
       FE(KKp)=AQ(KK1)
       QB(KKp)=AQ(KK2)
       IF (QB(KKp) .LE. 0.0) QB(KKp)=0.
       IF (FE(KKp) .LE. 0.0) FE(KKp)=0.
      ENDDO

      DO 600 KI=1,3
       KKI=(KI-1)*N
       KII=(KI-1)*KMAX1
       DO K=1,N
         QQ(K)=QQ1(KKI+K)
       ENDDO
       CALL COEFF (N,H,QQ,W,IOP,INT,WK)
       DO 63 K=1,KMAX1
        NB=0
        Y=Z(K)
       GO TO 66
   69  Y=Z1(K)
       QQ2(KII+K)=TAB(1)
   66  CALL TERP1 (N,H,QQ,W,Y,INT,TAB,ITAB)
       QQ3(KII+K)=TAB(1)
       NB=NB+1
       IF(NB.EQ.2) GO TO 63
       GO TO 69
   63 CONTINUE
  600 CONTINUE

      WRITE(6,1022) IVAA(2),(IVAA(I),I=5,7),IVAA(2),(IVAA(I),I=5,7)
      DO K=1,KMAX1
         WRITE(6,103) K,Z(K),(AA(K,I),I=1,3),Z1(K),(BB(K,I),I=1,3)
      ENDDO

      DO KK=1,KLES
         kkp=kk+1
         UB(KKp)=BB(KK,1)*100.
         VB(Kkp)=BB(KK,2)*100.
         WB(KKp)=AA(KK,3)
c        if (wb(kkp) .le. 0.0) wb(kkp)=0.
      ENDDO

      DO 700 KK=1,ITT
       do k=1,n
         dum(k)=uo_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 73 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        UB_6H(K+1,kk)=TAB(1)*100.
   73 CONTINUE
       do k=1,n
         dum(k)=vo_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 83 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        VB_6H(K+1,kk)=TAB(1)*100.
   83 CONTINUE
       do k=1,n
         dum(k)=wo_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 86 K=1,KLES
        Y=Z(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        WB_6H(K+1,kk)=TAB(1)
        if(y .ge. 16000.) WB_6H(K+1,kk)=0.
   86 CONTINUE
       do k=1,n
         dum(k)=q1o_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 88 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        Q1_6H(K+1,kk)=TAB(1)
        if(y .ge. 16000.) Q1_6H(K+1,kk)=0.
   88 CONTINUE
       do k=1,n
         dum(k)=q2o_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 89 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        Q2_6H(K+1,kk)=TAB(1)*.001
        if(y .ge. 16000.) q2_6H(K+1,kk)=0.
   89 CONTINUE


       do k=1,n
         dum(k)=t_adjust(k,kk)
       enddo
      call coeff(n,h,dum,w,iop,int,wk)
      do 173 k=1,kles
        y=z1(k)
        call terp1(n,h,dum,w,y,int,tab,itab)
        t_adjust0(k+1,kk)=tab(1)
  173 continue

       do k=1,n
         dum(k)=q_adjust(k,kk)
       enddo
      call coeff(n,h,dum,w,iop,int,wk)
      do 183 k=1,kles
        y=z1(k)
        call terp1(n,h,dum,w,y,int,tab,itab)
        q_adjust0(k+1,kk)=tab(1)
  183 continue


  700 CONTINUE

c      DO 75 K=1,KMAX
c         k1=k/2+1
c         k2=k/2
c        A1=AA(K1,1)
c         IF (K.EQ.K/2*2) A1=BB(K2,1)
c        A2=AA(K1,2)
c         IF (K.EQ.K/2*2) A2=BB(K2,2)
c       write(6,106) K,YYZ(K),AT(K),AQ(K),A1,A2
c   75 CONTINUE
      RETURN
  100 FORMAT(6X,I5,2F10.2,2x,F10.8,2x,F10.8)
  102 FORMAT(//,5X,6A8,A7,2A6)
 1021 FORMAT(//,6X,A8,A11,A10)
c1022 FORMAT(//,6X,A8,5A10)
 1022 FORMAT(//,6X,a8,5a10,a8,5A10)
  103 FORMAT(1X,I2,12F10.2)
 1031 FORMAT(I3,F9.2,F10.2,F8.2,F7.2,2F8.2,F7.2,2F6.2)
  106 FORMAT(2X,I4,8F10.3)
2031   format('                          in ainter',/,'          k     
     *z(k)     z1(k)     rm1(k)      rm(k)')
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TIMING(T)
      EXTERNAL ETIME         !workstation
      REAL*4 ETIME,TARRAY(2) !workstation
      T=ETIME(TARRAY)        !workstation	
      RETURN
      END
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine base (irs)
C     ******   SET DOMAIN AND INITIAL CONDITION
      PARAMETER (NX=130,NY=130,NZ=43,NM=NX,NT=2880,ITT=244)
      PARAMETER (NB1=10*NZ+5*NT)
      common/msui/ isounding,isui,kfine
      COMMON/OPTION/ LIPPS,IJKADV,ISTATMIN,IWATER,ITOGA,IMLIFTING,LIN,
     1   IRF,IADVH,IRFG,ISMG,ID

cccshie  11/16/01 add in for 3d!
      common/iceopt/ ice913,ilif

      common/option2/ isngtao,iwbar,iuvbar,isfc,ice,ice2,iradave,idelz
      COMMON/O4X/ A4H(NZ),A4K(NZ),D58X,D16X,D48X,D516X,D32X,D96X,D58Y,
     1   D16Y,D48Y,D516Y,D32Y,D96Y,D24X,D24Y
      COMMON/BXYZ/ IMAX,ILES,IL2,JMAX,JLES,JL2,KMAX,KLES,KL2,N,ISEC,
     1  NRAN,KT1,KT2
      COMMON/BX/ DX,D2X,RD2X,DX2,RDX,RD4X,RDX2,R2DX2,R4DX2
      COMMON/BY/ DY,D2Y,RD2Y,DY2,RDY,RD4Y,RDY2,R2DY2,R4DY2
      COMMON/BZ/ DZ,D2Z,RD2Z,DZ2,RDZ,RD4Z,RDZ2,R2DZ2,R4DZ2
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
      COMMON/B3UV/ U(NX,NY,NZ),V(NX,NY,NZ)
      COMMON/B4UV/ UU1(NX,NY,NZ),VV1(NX,NY,NZ)
      COMMON/DUMUW/ UMD(NX,NY,NZ),VMD(NX,NY,NZ),WMD(NX,NY,NZ)
      COMMON/B4/ TBSK(NZ),BSKT(NZ),BSKT2(NZ),BSKM(NZ),BSKM4(NZ),
     1   BSIT(NX,NY),BSIT2(NX,NY),BSIM(NX,NY),BSIM4(NX,NY),
     2   BSJT(NX,NY),BSJT2(NX,NY),BSJM(NX,NY),BSJM4(NX,NY)
      COMMON/B4Z/ BSITZ(NZ),BSIT2Z(NZ),BSIMZ(NZ),BSIM4Z(NZ),BSJTZ(NZ),
     1   BSJT2Z(NZ),BSJMZ(NZ),BSJM4Z(NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),TA(NZ),QA(NZ),TA1(NZ),
     1  QA1(NZ),COEF(NZ),C1(NZ),C2(NZ),C3(NZ),AM(NZ),AM1(NZ),UB(NZ),
     2  VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ),WBX(NX)
      COMMON/B6/ FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),ST(NZ),SV(NZ),
     1   SQ(NZ),SC(NZ),SE(NZ),SQA(NZ)
      common/bb6/ tls(nz),qls(nz),tls1(nz),tls2(nz),qls1(nz),qls2(nz),
     1   tls3(nz),tls4(nz),qls3(nz),qls4(nz),sft(nz),sfq(nz)  
      COMMON/CONT/ C38,C358,C610,C149,C879,C172,C409,C76,C218,C580,C141
      COMMON/BGWY/ UBAR,VBAR
      COMMON/DAMP/ RFA(NZ),RFA1(NZ),TBI(NZ),CNTC(NZ),CGWD
      COMMON/B555/ BA(NZ),BB(NZ)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/bb66/ ub_6h(nz,itt),ubt(nz),vb_6h(nz,itt),vbt(nz),
     1   wb_6h(nz,itt),wbt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     2   q2t(nz)
      common/q1q2z/ q1z_6h(nz,itt),q1tz(nz),q2z_6h(nz,itt),q2tz(nz)
      common/blsn/ ths(itt),qs(itt),ts(itt),
     1   pss(itt)
      COMMON/UPPER1/ T_ADJUST(NZ,ITT),Q_ADJUST(NZ,ITT),PRESS_T(NZ),
     1     TEMP_T(NZ),VAP_T(NZ)
      common/gbs/ tlsw(nz),qlsw(nz),ttlsw(nz),tqlsw(nz),smfff(nz),
     1   smffff(nz),tcof(nz),textra(nz),qextra(nz)
      common/bbb6/ wb1(nz),wb2(nz),wb3(nz)

      common/dinrad/ p00(nz),dz0(nz),tairsfc(nx,ny)
      common/dinradn/ qairsfc(nx,ny),pairsfc(nx,ny),thairsf(nx,ny)
      common/dinrad1/ dz1half
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON/BNDY2/ UB0(NZ),VB0(NZ)
      COMMON/BRH1/ SRR(NZ),QRR(NZ),Z11(NB1)
      COMMON/BSKD/ SKCOEF(NZ)
      COMMON/BCOR1/ DXR(NX,NY),DXXT(NX,NY),DXXT1(NX,NY),DYR(NX,NY),
     1   DYYT(NX,NY),DYYT1(NX,NY),DZR(NZ),DZZT(NZ),DZZT1(NZ)
      COMMON/BPBL/ UHT(NZ),WHT(NZ),TGBAT0
      COMMON/BA1/ Y1(NM),TAIR(NM),Z(NM),ZH(NM),QR(NM),QVS(NM),HU(NM),
     1   QVI(NM),HUI(NM),AN(NM),AN1(NM),AJ(NM),AJ1(NM),ANR(NM)
      COMMON/PICNST/ CPI

c      DIMENSION DZY(nz)
      REAL y2(nz)
c      DIMENSION DZYP(NZ)
      DIMENSION y3(nz),HGT1(NZ),HGT2(NZ),sst(itt),ssp(itt)
      DIMENSION z1(nz),z2(nz),y4(nz),hu1(nz)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      save
c      DATA DZY/0.,120.,150.,180.,210.,240.,270.,300.,350.,400.,450.,
c     1   500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,
c     2   1000.,1000.,1000.,1000.,1000.,1100.,1100.,1100.,1100.,1100.,
c     3   1200.,1200./
c      DATA DZYP/0.,120.,150.,180.,210.,240.,270.,300.,350.,400.,450.,
c     1   500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,
c     2   1000.,1000.,1000.,1000.,1000.,1100.,1100.,1100.,1100.,1100.,
c     3   1200.,1200./

      CPI=4.*ATAN(1.)
      ILES=NX-2
      JLES=NY-2
      KLES=NZ-2
      IL2=ILES
      JL2=JLES
      KL2=KLES
      ILES=ILES+1
      JLES=JLES+1
      KLES=KLES+1
      IMAX=ILES+1
      JMAX=JLES+1
      KMAX=KLES+1
      RIJL2=1./FLOAT(IL2*JL2)

      if(idelz.eq.1)then 
c        DZ=DZY(3)*100.
      else
        DZ1=20.      !m
      endif
      DX=2000.e2    !cm
      DY=2000.e2
      PSFC=1007.e3  
      TSFC=27.2     !C
      ktop100=29
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


      Z1(1)=0.
      Z2(1)=0.
       dze=dz1*100.
       do k=2,kmax
c        am(k)=1./(cc1+2.*cc2*(k-1.5)*dz1)
c        am1(k)=1./(cc1+2.*cc2*(k-2.)*dz1)
        y1(k)=(cc1+cc2*(k-1.5)*dz1)*(k-1.5)*dz1*100.
        y2(k)=(cc1+cc2*(k-2.)*dz1)*(k-2.)*dz1*100.
       enddo
        dz=dze
       do k=2,kmax
         z1(k)=y1(k)
         z2(k)=y2(k)
       enddo
      DO K=2,KLES
         AM(K)=DZ/(Z2(K+1)-Z2(K))
         AM1(K)=DZ/(Z1(K)-Z1(K-1))
      ENDDO
       am1(kmax)=dze/(y1(kmax)-y1(kles))
      AM1(1)=AM1(2)
      AM(1)=AM(2)
      AM(KMAX)=AM(KLES)
      do k=1,kmax
        hgt1(k)=z1(k)
        hgt2(k)=z2(k)
        uht(k)=.01*z1(k)
        wht(k)=.01*z2(k)
      enddo
      dz1half=z1(2)/100.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print*,'dz1half is',dz1half
c      write(6,*)' km      y11         y22        am        am1       
c     * z1          z2'
c      do k=kmax,1,-1
c       print*,z1(k)/100.,z2(k)/100.
c      enddo
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      D2X=DX+DX
      D2Y=DY+DY
      D2Z=DZ+DZ
      DX2=DX*DX
      DY2=DY*DY
      DZ2=DZ*DZ
      RDX=1./DX
      RDY=1./DY
      RDZ=1./DZ
      RD2X=.5*RDX
      RD2Y=.5*RDY
      RD2Z=.5*RDZ
      RD4X=.25*RDX
      RD4Y=.25*RDY
      RD4Z=.25*RDZ
      RDX2=RDX*RDX
      RDY2=RDY*RDY
      RDZ2=RDZ*RDZ
      R2DX2=.5*RDX2
      R2DY2=.5*RDY2
      R2DZ2=.5*RDZ2
      R4DX2=.25*RDX2
      R4DY2=.25*RDY2
      R4DZ2=.25*RDZ2
C     ****   4TH ORDER HORIZONTAL ADVECTION
      D16X=1./16.*RDX
      D48X=1./48.*RDX
      D58X=5./8.*RDX
      D32X=1./32.*RDX
      D96X=1./96.*RDX
      D24X=1./24.*RDX
      D516X=5./16.*RDX
      D16Y=1./16.*RDY
      D48Y=1./48.*RDY
      D58Y=5./8.*RDY
      D32Y=1./32.*RDY
      D96Y=1./96.*RDY
      D24Y=1./24.*RDY
      D516Y=5./16.*RDY
C     ******   RAYLEIGH RELAXTION TERMS ABOVE 15 KM
      RRTP0=10.
c     RRTP=16.
      RRTP=15.5  ! tao 5/1/01, cccshie 11/16/01 3d
      RRTPD=RRTP-RRTP0

      CGWD=1./(5.*24.*3600.)
      CG=1./(24.*3600.)
      CRF=1./900.

      DO K=1,KMAX
        A1=1.E-5*Z1(K)
        A2=1.E-5*Z2(K)
        RFA(K)=0.0
        RFA1(K)=0.0
        CNTC(K)=0.0
      IF(A1.GE.RRTP) CNTC(K)=CG
      IF(A1.LT.RRTP .AND. A1.GT.RRTP0) CNTC(K)=CG*(A1-RRTP0)/RRTPD
      IF(A1.GE.RRTP) RFA(K)=CRF*(A1-RRTP)
      IF(A2.GE.RRTP) RFA1(K)=CRF*(A2-RRTP)
      ENDDO
      CKH=2.
      CKHB=2.
      A2K=0.25E6
      AAK=0.001750/DT
c
      if (ijkadv .eq. 1) then
c        A2K=0.25E6
c        aak=0.00175/dt
         A2K=0.15E6 ! tao 5/1/01, cccshie 11/16/01 3d
         aak=0.001250/dt  ! tao 5/1/01, cccshie 11/16/01 3d
      endif
c      kref=7
C      if (kfine .eq. 0) kref=2
C      if (idelz .eq. 1) kref=18
      zkref=250.    !m
      kref=2
      do k=2,kmax
        if(((z2(k)-z2(k-1))/100.).le.zkref) kref=k
      enddo
      print*,'kref is',kref
      A2I=AAK*DX*DX
      A2J=AAK*DY*DY
      DO 20 K=1,KMAX
        krefh=kref
        Y3(K)=DZ*DZ/(AM(K)*AM(K))
        COEF(K)=CK*CK*Y3(K)*.5
        C1(K)=3.*980.*COEF(K)
        C2(K)=.5*CE/(CK*Y3(K))
        C3(K)=2./(3.*CK*CK*Y3(K))
ctao  (11-13)
        A4H(K)=AAK*(AM(kref)/AM(K))
c        A4H(K)=AAK
          A4K(K)=A4H(K)/CKHB
        TBSK(K)=2.0e6*(AM(kref)/AM(K)) ! tao 5/1/01, cccshie 11/16 3d
c       TBSK(K)=3.0e6*AM(kref)/AM(K)
             SKCOEF(K)=0.025*DX*DX
        BSKT(K)=A2K*AM(kref)/AM(K)
          BSKT2(K)=2.*BSKT(K)
        BSKM(K)=BSKT(K)/CKH
          BSKM4(K)=4.*BSKM(K)
c
        BSITZ(K)=A2I*AM(krefh)/AM(K)
          BSIT2Z(K)=2.*BSITZ(K)
        BSIMZ(K)=BSITZ(K)/CKH
          BSIM4Z(K)=4.*BSITZ(K)/CKH
        BSJTZ(K)=A2J*AM(krefh)/AM(K)
          BSJT2Z(K)=2.*BSJTZ(K)
        BSJMZ(K)=BSJTZ(K)/CKH
          BSJM4Z(K)=4.*BSJTZ(K)/CKH
   20 CONTINUE
      DO J=1,JMAX
      DO I=1,IMAX
       BSIT(I,J)=A2I
       BSIT2(I,J)=2.*A2I
       BSIM(I,J)=A2I/CKHB
       BSIM4(I,J)=4.*A2I/CKHB
       BSJT(I,J)=A2J
       BSJT2(I,J)=2.*A2J
       BSJM(I,J)=A2J/CKHB
       BSJM4(I,J)=4.*A2J/CKHB
      ENDDO
      ENDDO
      if(isounding.eq.1) then 
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
c        call ainter(sst,ssp,tb,fd,qb,fe,ub,vb,wb,z1,z2,dz,topht)
        call ainter(sst,ssp,tb,fd,qb,fe,ub,vb,wb,z1,z2,topht)
        PSFC=SSP(1)*1000.
         do k=1,kmax
           if((z2(k)/100.).lt.topht) ktop100=k+1
cccshie 11/16/01
          if((z2(k)/100.).le.topht) ktoptq=k
          if((z2(k)/100.).le.topht) ktopuv=k

         enddo
       print*,'ktop100 is',ktop100,.01*z2(ktop100),.01*z1(ktop100),topht

cccshie 11/16/01
       print*,'ktoptq, m, topht ',ktoptq,.01*z2(ktoptq),topht ! toptq use topht
       print*,'ktopuv, m, topht ',ktopuv,.01*z2(ktopuv),topht ! topuv use topht

       endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,210)
      DO K1=2,KMAX
        K=KMAX+2-K1
c     WRITE(6,2031) K,TB(K),QB(K),FD(K),FE(K),TLS(K),QLS(K),QR(K),
      WRITE(6,2031) K,TB(K),QB(K),FD(K),FE(K),QR(K),  ! cccshie 11/16/01
     1               UB(K),VB(K),WB(K)
      ENDDO
      WRITE(6,204)
      DO K1=2,KMAX
        K=KMAX+2-K1
      WRITE(6,2441) K,Z1(K)/100.,Z2(K)/100.,AM(K),AM1(K),COEF(K)
     *,              C1(K),C2(K),C3(K),TBSK(K),BSKM(K),RFA(K) ! 11/16/01 TBSK,BSKM
      ENDDO
C     *********************
      DO K=2,KMAX
        FD(K)=FD(K)+273.16
        FE(K)=FE(K)*1.E-3
        TB(K)=TB(K)+273.16
        QB(K)=QB(K)*1.E-3
      ENDDO
      TGBAT0=FD(2)+0.75           !sea surface temperature 0.75 degree warmer

      A1=980.*DZ/RA
      P0(2)=LOG(PSFC)
      DO K=3,KMAX
        P0(K)=P0(K-1)-A1/(TB(K-1)*(1.+.61*QB(K-1))*AM(K-1))
      ENDDO
      RHO1(KMAX)=EXP(P0(KMAX))/(RA*FD(KMAX)*(1.+.61*FE(KMAX)))
      DO  K=2,KLES
        DZ0(K)=Z2(K+1)-Z2(K)
      ENDDO
      DZ0(KMAX)=DZ0(KLES)
      DO K=2,KMAX
         P00(K)=EXP(P0(K))
      ENDDO
      DO 50 K=2,KLES
        RHO1(K)=EXP(P0(K))/(RA*FD(K)*(1.+.61*FE(K)))
        P0(K)=EXP((P0(K)+P0(K+1))*.5)
        PI(K)=(P0(K)/1000.E3)**.286
        F0(K)=AL/(PI(K)*CP)
        Y1(K)=TB(K)
        TB(K)=TB(K)/PI(K)
        QVS(K)=3.799E3/P0(K)*EXP(17.26939*(Y1(K)-273.16)/(Y1(K)-35.86))
        Y3(K)=3.799E3/P0(K)*EXP(21.87456-5807.695/(Y1(K)-7.66))

cccshie 11/16/01 add in for 3d
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
      RHO(K)=P0(K)/(RA*Y1(K)*(1.+.61*QB(K)))
   50 CONTINUE
      DO K=3,KLES-1
        IF(TB(K).LT.TB(K-1))TB(K)=.5*(TB(K-1)+TB(K+1))
        IF(TB(K).LT.TB(K-1)) TB(K)=TB(K-1)-TB(K)+TB(K-1) ! shie 11/16/01 3d
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
            q1_6h(k,icy)=q1_6h(k350,icy)*divfac  ! 11/26/01 shie found typo!
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
c       if(iuvbar.eq.0)then
c         A2=0.
c         A3=0.
c         A4=0.
c         DO 6666 K=2,KLES
c           A1=RHO(K)*DZ/AM(K)
c           A2=A2+A1*UB(K)
c           A3=A3+A1*VB(K)
c           A4=A4+A1
c 6666   CONTINUE
c         DO 6669 K=2,KLES
c 6669     VB(K)=A3/A4
c       endif

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
          THS(k)=TSL
          QS(k)=QSL
          TS(k)=tsfc+273.16
          PSS(k)=psfc
        enddo
      enddo

        write(6,12121) k,tsl,qsl*1000.,psfc/1000.,tsfc+273.16
12121 format(2x,'it=',i4,4x,'tsfc=',f12.5,4x,'qsfc=',f12.5,4x,
     1          'psfc=',f12.5,4x,'tsfc c=',f12.5)
      TB(1)=TB(2)
      QB(1)=QB(2)
      UB(1)=UB(2)
      VB(1)=VB(2)
      WB(1)=0.
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
        if(k.le.2)wb(k)=0.
        if (k .eq. kmax) wb(k)=0.
          if (p0(k) .le. 75.e3) then
            if (wb(k) .lt. 0.0) wb(k)=0.
          endif
        do icy=1,itt
          wb_6h(k,icy)=scale*wb_6h(k,icy)/(rho1(k)*980.)
          if (k .le. 2) wb_6h(k,icy)=0.
          if (k .eq. kmax) wb_6h(k,icy)=0.
            if (p0(k) .le. 75.e3) then
              if (wb_6h(k,icy) .lt. 0.0) wb_6h(k,icy)=0.
            endif
        end do
        RRHO(K)=1./RHO(K)
        RRHO1(K)=1./RHO1(K)
        TBI(K)=TB(K)
        UB(K)=UB(K)
        VB(K)=VB(K)
        RRHO(K)=1./RHO(K)
        RRHO1(K)=1./RHO1(K)
        SRR(K)=1./SQRT(RHO(K))
        QRR(K)=SQRT(SRR(K))
        UB0(K)=UB(K)
        VB0(K)=VB(K)
   70 CONTINUE
c$doacross local(j,i)
       DO J=1,JMAX
       DO I=1,IMAX
         DXR(I,J)=RDX
         DXXT1(I,J)=DT*RDX
         DXXT(I,J)=DT*RDX
       ENDDO
       ENDDO
       DO K=2,KLES
        DZR(K)=AM(K)*RDZ/RHO(K)
        DZZT1(K)=DT*RDZ*AM(K)
        DZZT(K)=DT*RDZ*AM1(K)
       ENDDO
        DZR(1)=DZR(2)
        DZZT(1)=DZZT(2)
        DZR(KMAX)=DZR(KLES)
        DZZT(KMAX)=DZZT(KLES)
c$doacross local(j,i)
       DO J=1,JMAX
       DO I=1,IMAX
        DYR(I,J)=RDY
        DYYT1(I,J)=DT*RDY
        DYYT(I,J)=DT*RDY
       ENDDO
       ENDDO
      IF(IRS.EQ.1) GO TO 1111
      DO K=1,KMAX
        TA(K)=TB(K)
        QA(K)=QB(K)
        TA1(K)=TB(K)
        QA1(K)=QB(K)
        UB1(K)=UB(K)
        VB1(K)=VB(K)
c$doacross local(j,i)
      DO J=1,JMAX
      DO I=1,IMAX
        UU1(I,J,K)=UB(K)
        VV1(I,J,K)=VB(K)
        U(I,J,K)=UB(K)
        V(I,J,K)=VB(K)
        umd(I,J,K)=UB(K)
        vmd(i,J,k)=vb(k)
        wmd(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
 1111 WRITE(6,200)
      DO K1=2,KLES
      K=KLES+2-K1
      KM=K-1
      WRITE(6,2032)KM,P0(K)/1000.,PI(K),TB(K),QB(K)*1000.,QVS(K)*1000.
     *,     HU1(K),BA(K)*1000.,BB(K)*1000.,Y3(K)*1000.,Y4(K)
      ENDDO
      WRITE(6,201)
      DO K1=2,KLES
      K=KLES+2-K1
      KM=K-1
c     WRITE(6,2033)KM,RHO(K)*1000.,RHO1(K)*1000.,TLS(K),QLS(K),F0(K)
      WRITE(6,2043)KM,RHO(K)*1000.,RHO1(K)*1000.,F0(K) ! cccshie 11/16/01 for 3d
     *,FD(K),FE(K),UB(K)/100.,VB(K)/100.,WB(K),BSIMZ(K)
      ENDDO
      do icy=1,itt
        print*,'          '
        print*,'itt =',icy
        print*,'          '
        write(6,241)
        do k1=2,kles
          k=kles+2-k1
          km=k-1
          write(6,203) km,P0(K)/1000.,ub_6h(k,icy),vb_6h(k,icy),
     1                   wb_6h(k,icy),Q1_6H(k,icy),Q2_6H(k,icy)*1000.,
     2                   t_adjust(k,icy),q_adjust(k,icy)*1000.
        enddo
      enddo
      RETURN
  200 FORMAT(//,' LEVEL   PRES    PI      TB       QB     QSW     HUM      
     *   BA      BB     QSI      HUI')
c 201 FORMAT(//,'LEVEL  RHO   RHO1    TLS     QLS     F0        FD     
c    *    FE       UB      VB    WB      BSIM')
cccshie 11/16/01 for 3d
c 201	FORMAT(//,'LEVEL  RHO   RHO1      F0        FD          FE       UB      VB    WB      BSIM')
  201	FORMAT(//,'LEVEL  RHO   RHO1      F0        FD          FE       
     1UB      VB    WB      BSIM')
  202 FORMAT(16F5.2)
  203 FORMAT(I4,10F12.5,/)
 2035 FORMAT(I4,2F10.5,4F12.3)
 2031 FORMAT(1X,I4,2(F8.2,F7.2),2F8.2,F7.2,2F10.2,F8.2)
 2032 FORMAT(I4,F9.2,F7.3,F9.2,7F8.2)
 2033 FORMAT(I4,2F7.3,F8.3,F7.3,F9.2,E11.3,E12.3,2F7.2,F7.3,E11.3)
 2043 FORMAT(I4,2F7.3,F9.2,E11.3,E12.3,2F7.2,F7.3,E11.3) ! cccshie 11/16/01
  204 FORMAT(//,'LEVEL  Z      Z1     AM    AM1    COEF       C1       
     *C2       C3        TBSK       BSKM')
  241 format(//,'level          p           u           v           w 
     1        q1          q2       tdt       qdt')
 2441 FORMAT(I3,2F8.1,2F6.3,5E10.3,E9.2,E10.3)
  210 FORMAT(//,'  LEVEL  TB      QB     FD     FE1      QR       UB
     1     VB        WB')
  250 format(16f5.0)
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine obs_data(sst,ssp,press,temp,vap,uu,vv,ww
     *,                  ubo_6h,vbo_6h,wbo_6h,q1o_6h,q2o_6h)
c     parameter(ITT=244,npin=33,ITTSKIP=16) ! cccshie 11/15/01 cause problem for scsmex
      parameter(ITT=244,npin=41,ITTSKIP=0)
      COMMON/BT/ DT,D2T,RIJL2,DTS,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,
     1    PSFC,FCOR,SEC,AMINUT,RDT
      dimension press(npin),temp(npin),vap(npin),
     *     uu(npin),vv(npin),ww(npin)
     *,    ubo_6h(npin,itt),vbo_6h(npin,itt),wbo_6h(npin,itt)
     *,    q1o_6h(npin,itt),q2o_6h(npin,itt)
      common/upper/ t_adjust(npin,itt),q_adjust(npin,itt),press_t(npin),
     1      temp_t(npin),vap_t(npin)
      dimension ttt(npin,itt),qqq(npin,itt),ppp(npin,itt),sst(itt),
     1  ssp(itt)

      real dum(npin),q1horz(npin),q2horz(npin),q1vert(npin),q2vert(npin)
      real h(npin)

      save

c     open(60,file='/usr/raid1/djohnson/case3/ifa_sst.dat'
c     1       ,status='old')
c      open(61,file='/usr/raid1/djohnson/case3/ifa_sounding.dat'
c     1       ,status='old')
c      open(62,file='/usr/raid1/djohnson/case3/ifa_surface.dat'
c     1       ,status='old')
c      open(63,file='/usr/raid1/djohnson/case3/ifa_forcing.dat'
c     1       ,status='old')

      open(60,file='input',status='old')
      open(61,file='t161q161.out',status='old')

C  Read into arrays, the times at beginning of simulation to be written over


        READ(60,202) (dum(K),dum(K),K=2,NPIN)
        READ(60,202) (dum(K),dum(K),K=2,NPIN)
        READ(60,250) (dum(K),K=2,NPIN)
        READ(60,250) (dum(K),K=2,NPIN)
        READ(60,202) (dum(K),K=2,NPIN)

          DO ICY=1,ittskip
            READ(60,250) (dum(K),K=2,NPIN)
            READ(60,250) (dum(K),K=2,NPIN)
            READ(60,202) (dum(K),K=2,NPIN)
          ENDDO


C  Now read in the data we want to use in the simulation


        do i=1,itt
         READ(60,250) (ubo_6H(K,I),K=2,NPIN)
         READ(60,250) (vbo_6H(K,I),K=2,NPIN)
         READ(60,202) (wbo_6H(K,I),K=2,NPIN)
         do k=2,npin
            ubo_6H(K,I)=ubo_6h(k,i)*.01
            vbo_6H(K,I)=vbo_6h(k,i)*.01
         enddo
        enddo

      rewind(60)
      read(60,202) (dum(k),dum(k), k=2,npin)
      read(60,202) (dum(k),dum(k), k=2,npin)
      read(60,250) (dum(k), k=2,npin)
      read(60,250) (dum(k), k=2,npin)
      read(60,202) (dum(k), k=2,npin)

      do i=1,161
          read(60,250) (dum(k), k=2,npin)
          read(60,250) (dum(k), k=2,npin)
          read(60,250) (dum(k), k=2,npin)
      enddo

          DO I=1,ittskip
          read(60,202) (dum(k), k=2,npin)
          read(60,202) (dum(k), k=2,npin)
          read(60,202) (dum(k), k=2,npin)
          read(60,202) (dum(k), k=2,npin)
          ENDDO

          DO I=1,itt
          read(60,202) (q1horz(k), k=2,npin)
          read(60,202) (q1vert(k), k=2,npin)
          read(60,202) (q2horz(k), k=2,npin)
          read(60,202) (q2vert(k), k=2,npin)

          do k=2,npin
             q1o_6h(k,i)=q1horz(k)+q1vert(k)
             q2o_6h(k,i)=(-q2horz(k)-q2vert(k))/2.49
          enddo
          ENDDO

        do i=1,ittskip
           READ(61,202) (dum(K),K=2,npin)
           READ(61,202) (dum(K),K=2,npin)
        enddo

          DO I=1,itt
           READ(61,202) (ttt(k,i),K=2,npin)
           READ(61,202) (qqq(k,i),K=2,npin)
          ENDDO


      do i=1,itt
        ppp(1,i)=1007.
        cc1=10.
        cc2=.01
        dz=34.e2
        H(1)=0.
         DO K=2,NPIN
             H(k)=(cc1+cc2*(k-1.5)*dz*.01)*(k-1.5)*dz
             h(k)=h(k)*.01
c             Z2(k)=(cc1+cc2*(k-2.)*dz*.01)*(k-2.)*dz
         ENDDO

         H(1)=0.
         TERM1=2.87E+6/980.*.5

         DO K=2,NPIN
            KM=K-1

       ppp(k,i)=ppp(km,i)
     1   *exp((h(km)-h(k))*100./TERM1/((ttt(K,i)+273.16)
     1   *(1.+.61*qqq(K,i)*1.E-3)+(ttt(KM,i)+273.16)
     1   *(1.+.61*qqq(KM,i)*1.E-3)))

        enddo
        enddo


      do 999 i=1,itt

        PRINT*,'NON-MODIFIED LARGE SCALE CONDITIONS AT',I*3-3,'HOURS'
        WRITE(6,240)
        DO K1=1,NPIN-1
          K=NPIN+1-K1
          WRITE(6,200) K,PPP(K,I),TTT(K,I),QQQ(K,I),UBO_6H(K,I),
     2                 VBO_6H(K,I),WBO_6H(K,I),Q1O_6H(K,I),Q2O_6H(K,I)
        ENDDO

c       READ(60,*) SST(I)
c       READ(62,*) SSP(I)

        sst(i)=27.2
        ssp(i)=1007.

        PPP(1,I)=SSP(1)
c       PPP(1,I)=1007.
c       PPP(NPIN-1,I)=50.
c       PPP(NPIN,I)=25.
        PSFC=SSP(1)*1000.

        TTT(1,I)=((TTT(2,I)+273.16)*(PPP(1,I)/PPP(2,I))**(RA/CP))-273.16
        IF (TTT(NPIN-2,I)-TTT(NPIN-3,I) .GT. 10.)
     1                                   TTT(NPIN-2,I)=TTT(NPIN-3,I)+10.
c       TTT(NPIN-1,I)=TTT(NPIN-2,I)+10.
c        TTT(NPIN,I)=TTT(NPIN-1,I)+10.

        QQQ(1,I)=QQQ(2,I)
        QQQ(NPIN-1,I)=0.
        QQQ(NPIN,I)=0.

c        DO K=35,NPIN-3
c          IF (UBO_6H(K,I) .LE. -20.0) UBO_6H(K,I)=-20.0
c        ENDDO

c        DO K=2,NPIN-2
c          IF (PPP(K,I) .Le. 100.) Q1O_6H(K,I)=0.
c          IF (PPP(K,I) .Le. 100.) Q2O_6H(K,I)=0.
c         IF (PPP(K,I) .Le. 100.) QQQ(K,I)=0.
c          IF (PPP(K,I) .Le. 100.) WBO_6H(K,I)=0.0
c          if (i .ge. 5) then
c            IF (PPP(K,I) .LE. 125.) THEN
c             WBO_6H(K,I)=0.75*((REAL(NPIN)-REAL(K))/9.)**2*WBO_6H(K,I)
c             Q1O_6H(K,I)=0.5*Q1O_6H(K,I)
c             Q2O_6H(K,I)=0.5*Q2O_6H(K,I)
c             if (Q1O_6H(K,I) .ge. 1.5) Q1O_6H(K,I)=1.5
c             if (Q1O_6H(K,I) .le. -1.5) Q1O_6H(K,I)=-1.5
c            ENDIF
c          else
c            IF (PPP(K,I) .LE. 150.) THEN
c             WBO_6H(K,I)=0.
c             Q1O_6H(K,I)=0.
c             Q2O_6H(K,I)=0.
c            ENDIF
c          endif
c        ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        dpp=(ppp(1,i)-ppp(2,i))/(ppp(2,i)-ppp(3,i))
c        ubo_6h(1,i)=ubo_6h(2,i)+(ubo_6h(2,i)-ubo_6h(3,i))*dpp
c        if ( abs(ubo_6h(1,i)) .ge. abs(ubo_6h(2,i)) )
c     1      ubo_6h(1,i)=ubo_6h(2,i)
c        ubo_6h(npin-3,i)=ubo_6h(npin-4,i)
c        ubo_6h(npin-2,i)=ubo_6h(npin-3,i)
c        ubo_6h(npin-1,i)=ubo_6h(npin-2,i)
c        ubo_6h(npin,i)=ubo_6h(npin-1,i)

c        vbo_6h(1,i)=vbo_6h(2,i)+(vbo_6h(2,i)-vbo_6h(3,i))*dpp
c        if ( abs(vbo_6h(1,i)) .ge. abs(vbo_6h(2,i)) )
c     1     vbo_6h(1,i)=vbo_6h(2,i)
c        vbo_6h(npin-3,i)=vbo_6h(npin-4,i)
c        vbo_6h(npin-2,i)=vbo_6h(npin-3,i)
c        vbo_6h(npin-1,i)=vbo_6h(npin-2,i)
c        vbo_6h(npin,i)=vbo_6h(npin-1,i)

        wbo_6h(1,i)=0.
        wbo_6h(npin-1,i)=0.
        wbo_6h(npin,i)=0.

c        q1o_6h(1,i)=q1o_6h(2,i)+(q1o_6h(2,i)-q1o_6h(3,i))*dpp
c        if ( abs(q1o_6h(1,i)) .ge. abs(q1o_6h(2,i)) )
c     1     q1o_6h(1,i)=q1o_6h(2,i)
c        q1o_6h(npin-1,i)=0.
c        q1o_6h(npin,i)=0.

c        q2o_6h(1,i)=q2o_6h(2,i)+(q2o_6h(2,i)-q2o_6h(3,i))*dpp
c        if ( abs(q2o_6h(1,i)) .ge. abs(q2o_6h(2,i)) )
c     1     q2o_6h(1,i)=q2o_6h(2,i)
c        q2o_6h(npin-1,i)=0.
c        q2o_6h(npin,i)=0.

999   continue

      DO K=1,NPIN
        PRESS_T(K)=0.
        TEMP_T(K)=0.
        VAP_T(K)=0.
      ENDDO
      SSP_T=0.
      SST_T=0.

      DO I=1,ITT
        IP=I+1
        IF(I.EQ.ITT) IP=ITT
        SSP_T=SSP_T+SSP(I)
        SST_T=SST_T+SST(I)
        DO K=1,NPIN
          PRESS_T(K)=PRESS_T(K)+PPP(K,I)
          TEMP_T(K)=TEMP_T(K)+TTT(K,I)
          VAP_T(K)=VAP_T(K)+QQQ(K,I)
          T_ADJUST(K,I)=TTT(K,IP)-TTT(K,I)
          Q_ADJUST(K,I)=QQQ(K,IP)-QQQ(K,I)
        ENDDO
      ENDDO

      DO K=1,NPIN
        PRESS_T(K)=PRESS_T(K)/ITT
        TEMP_T(K)=TEMP_T(K)/ITT
        VAP_T(K)=VAP_T(K)/ITT
      ENDDO
      SSP_T=SSP_T/ITT
      SST_T=SST_T/ITT

      PRINT*
      PRINT*,'MEAN SOUNDING'
      PRINT*
      WRITE(6,241)
      DO K1=1,NPIN
        K=NPIN+1-K1
        WRITE(6,200) K,PRESS_T(K),TEMP_T(K),VAP_T(K)
      ENDDO

      DO K=1,NPIN
        PRESS(K)=PPP(K,1)
        TEMP(K)=TTT(K,1)
        VAP(K)=QQQ(K,1)
        UU(K)=UBO_6H(K,1)
        VV(K)=VBO_6H(K,1)
        WW(K)=WBO_6H(K,1)
      ENDDO
      print*
      print*,'initial sounding'
      print*
      write(6,241)
      do k1=1,npin
        k=npin+1-k1
       write(6,200) k,press(k),temp(k),vap(k),uu(k),vv(k),ww(k)
      enddo

      do icy=1,itt
        print*
        print*,'itt =',icy
        WRITE(6,2404)
        do k1=1,npin
          K=npin+1-K1
       WRITE(6,200) K,ppp(K,icy),ttt(K,icy),qqq(K,icy),T_ADJUST(K,Icy),
     1                 Q_ADJUST(K,Icy),UBo_6H(K,ICY),VBo_6H(K,ICY),
     2                 WBo_6H(K,ICY),Q1o_6H(K,ICY),Q2o_6H(K,ICY)
        enddo
      enddo
  202 FORMAT(16F5.2)
  250 format(16f5.0)

      return

  200 format(i4,10f12.5)
  240 format(//,'level          p           t           q         u
     1     v            w           q1        q2')
 2404 format(//,'level          p           t           q         tdt
     1    qdt          u            v          w          q1        q2')
  241 format(//,'level          p           t           q           u
     1       v           w')
      end

c-----------------------------------------------------------------------
      subroutine obs(sst,ssp,press,temp,vap,uu,vv,ww,ubo_6h,vbo_6h
     1                    ,wbo_6h,q1o_6h,q2o_6h)
      parameter(NPIN=41,ITT=244,ITTSKIP=0)

      dimension ttt(npin,itt),qqq(npin,itt),ppp(npin,itt),sst(itt),
     1          ssp(itt),theta_6h(npin,itt),div(npin,itt)

      common /zzzobs/ zzz(npin,itt)

      dimension hu(npin,itt),vu(npin,itt),hv(npin,itt),vva(npin,itt),
     1          ht(npin,itt),vt(npin,itt),hq(npin,itt),vq(npin,itt)

      dimension press(npin),temp(npin),vap(npin),uu(npin),vv(npin),
     1          ww(npin),ubo_6h(npin,itt),vbo_6h(npin,itt),
     2          wbo_6h(npin,itt),q1o_6h(npin,itt),q2o_6h(npin,itt)
      common/upper/ t_adjust(npin,itt),q_adjust(npin,itt),press_t(npin),
     1      temp_t(npin),vap_t(npin)
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      parameter(ISKIPMAY=0)
cccshie3
c     parameter(ISTARTDAY=6)  ! 4/20 rerun,3/16/01 shie, skip 5 days for 1st 10days
c     parameter(ISTARTDAY=15) ! 4/18/01,4/6/01 shie, skip 17 days for 1st 10days

      parameter(ISTARTDAY=18)  ! 10/23/01, 051800-052600 ( 0+18=18)
c     parameter(ISTARTDAY=33)  ! 10/23/01, 060200-061100 (31+ 2=33)

      parameter(ISTARTHOUR=0) ! 4/20, 4/18, 4/8

cccshie average over the troble data at 052018(may 20, 18z):  iavge=1
c      iavge=1  ! average over the trouble data at 052018 (old data only)
       iavge=0  ! use the original trouble data at 052018, or use new data 6/5/01

      open(60,file='fields.nesa',status='old')
      open(61,file='advect.nesa',status='old')

C  Read into arrays, the times at beginning of simulation to be written over

      iskip = (iskipmay*31 + istartday-1) * 4  + istarthour/6

      if(iskip.gt.0)then
         do i=1,iskip
            read(60,*)
            read(61,*)
	    do k=1,npin
	       read(60,*) 
	       read(61,*) 
	    enddo
         enddo
      endif

C  Now read in the data we want to use in the simulation


      do 999 iii=iskip+1,itt
         i=iii-iskip
         read(60,*)
         read(61,*)

	 do k=1,npin

	  read(60,*) ppp(k,i),zzz(k,i),ubo_6h(k,i),vbo_6h(k,i),
     1               wbo_6h(k,i),ttt(k,i),theta_6h(k,i),qqq(k,i),
     1               div(k,i)
	  read(61,*) ppp(k,i),hu(k,i),vu(k,i),hv(k,i),vva(k,i),ht(k,i),
     1               vt(k,i),hq(k,i),vq(k,i)

ccc  03-09-2001 Change, dan and shie

cccshie !!!! Important!!!
c
c 3/30/01 shie summarize (hope this should be the correct arrangement after all)
c
c  since the large-scale advection terms will be substracted later, so the terms
c  dealt here should be changed to negative of them
c
c   (1) the large-scale Q advections 
c    "-(W(dQ/dz)+U(dQ/dx))" = -(vq(k,i)+hq(k,i))
c     so (chnage sign) -->
c    "W(dQ/dz)+U(dQ/dx)" = (vq(k,i)+hq(k,i)), which are assigned to "q2o_6h(k,i)"
c
c   (2) the large-scale T advections 
c    "-(W(dT/dz)+U(dT/dx))" = "-omega(dT/dp)+omega/(rho*cp)-U(dT/dx)"
c     = -vt(k,i)+omega/(rho*cp)-ht(k,i))
c     so (chnage sign) -->
c    "W(dT/dz)+U(dT/dx)" = "omega(dT/dp)-omega/(rho*cp)+U(dT/dx)"
c     = (vt(k,i)-omega/(rho*cp)+ht(k,i)), which are assigned to "q1o_6h(k,i)"
c
cccshie !!!! Important!!!

          density=ppp(k,i)*100./(287.*(ttt(k,i)+273.16)
     1           *(1.+.61*qqq(k,i)*.001))                   ! S.I units
c
c "vt(k,i)-omega/(rho*cp)" and convert units (mb-->pascal; hr-->sec)
c
          vert_for=vt(k,i)-100./3600.*wbo_6h(k,i)/1004./density  ! omega/cp/rho

          q1o_6h(k,i)=(ht(k,i)+vert_for)*24.*3600.   ! correct, 3/30/01
c            print*,'q1o',i,k,q1o_6h(k,i),wbo_6h(k,i)
c         q1o_6h(k,i)=-(ht(k,i)+vert_for)*24.*3600.  ! wrong, 3/30/01,temp LSF
cc	  q1o_6h(k,i)=(ht(k,i)+vt(k,i))*24.*3600.    ! wrong, 3/13/01 summarize

ccc
c	  q2o_6h(k,i)=-(hq(k,i)+vq(k,i))*24.*3600. ! wrong 3/30/01 moisture LSF
	  q2o_6h(k,i)=(hq(k,i)+vq(k,i))*24.*3600.  ! correct 3/30/01 moisture LSF
	enddo

c	PRINT*,'NON-MODIFIED LARGE SCALE CONDITIONS AT',I*6-6,'HOURS'
c	WRITE(6,210)
c	DO K1=1,NPIN-1
c	  K=NPIN+1-K1
c	  WRITE(6,200) K,PPP(K,I),TTT(K,I),QQQ(K,I),UBO_6H(K,I),
c     2                 VBO_6H(K,I),WBO_6H(K,I),Q1O_6H(K,I),Q2O_6H(K,I)
c	ENDDO

	ssp(i)=ppp(1,i)
	sst(i)=ttt(1,i)


999   continue

cccshie average over the trouble data at 052018(may 20, 18z):  iavge=1
       if(iavge.eq.1) then
       do k=1,npin
          q1o_6h(k,12)=(q1o_6h(k,11)+q1o_6h(k,13))*.5
          wbo_6h(k,12)=(wbo_6h(k,11)+wbo_6h(k,13))*.5
c            print*,'q1o1',i,k,q1o_6h(k,12),wbo_6h(k,12)
       enddo 
       endif

      TSFC=SST(1)
      PSFC=SSP(1)*1000.

      DO K=1,NPIN
        PRESS_T(K)=0.
        TEMP_T(K)=0.
        VAP_T(K)=0.
      ENDDO

      do iii=iskip+1,itt
         i=iii-iskip
        DO K=1,NPIN
          PRESS_T(K)=PRESS_T(K)+PPP(K,I)
          TEMP_T(K)=TEMP_T(K)+TTT(K,I)
          VAP_T(K)=VAP_T(K)+QQQ(K,I)
        ENDDO
      ENDDO


      DO K=1,NPIN
        PRESS_T(K)=PRESS_T(K)/(ITT-ISKIP)  ! cccshie note "iskip" added by dan
        TEMP_T(K)=TEMP_T(K)/(ITT-ISKIP)
        VAP_T(K)=VAP_T(K)/(ITT-ISKIP)
      ENDDO


      PRINT*
      PRINT*,'MEAN SOUNDING'
      PRINT*
      WRITE(6,241)
      DO K1=1,NPIN
        K=NPIN+1-K1
        WRITE(6,200) K,PRESS_T(K),TEMP_T(K),VAP_T(K)
      ENDDO


      DO K=1,NPIN
        PRESS(K)=PPP(K,1)
        TEMP(K)=TTT(K,1)
        VAP(K)=QQQ(K,1)
        UU(K)=UBO_6H(K,1)
        VV(K)=VBO_6H(K,1)
        WW(K)=WBO_6H(K,1)
      ENDDO

      print*
      print*,'initial sounding'
      print*
      write(6,241)
      do k1=1,npin
        k=npin+1-k1
       write(6,200) k,press(k),temp(k),vap(k),uu(k),vv(k),ww(k)
      enddo

      do icy=1,itt
        print*
        print*,'itt =',icy
        WRITE(6,240)
        do k1=1,npin
          K=npin+1-K1
       WRITE(6,200) K,ppp(K,icy),ttt(K,icy),qqq(K,icy),T_ADJUST(K,Icy),
     1                 Q_ADJUST(K,Icy),UBo_6H(K,ICY),VBo_6H(K,ICY),
     2                 WBo_6H(K,ICY),Q1o_6H(K,ICY),Q2o_6H(K,ICY)
        enddo
      enddo

      return

  200 format(i4,10f12.5)
  210 format(//,'level        p           t           q          u   
     1      v           w         q1        q2')
  240 format(//,'level          p           t           q         tdt 
     1     qdt           u           v           w         q1        q2
     2 ')
  241 format(//,'level          p           t           q           u 
     1       v           w')
      end
c-----------------------------------------------------------------------
