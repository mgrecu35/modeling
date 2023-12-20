      subroutine consatl (rho)
C     (LIN) SPECIFY SOME CONSTANTS IN SATICE ROUTINE   ******
      parameter (NZ=43,NT=2880)
      parameter (nb=10*nz+5*nt)

      real    rho(nz)
      real    srr(nz),qrr(nz),z1(nb)
      common/brh1/ srr,qrr,z1

      common/ilin/ lin

      common/iceopt/ ice913,ilif

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141

      real    ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq

      real    tnw,tns,tng,roqs,roqg,roqr
      common/size/ tnw,tns,tng,roqs,roqg,roqr

      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC
      COMMON/BSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,BND1,RN2,BND2,
     1   RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,RN101,RN10A,RN11,RN11A,
     2   RN12,RN12A(31),RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN17,
     3   RN17A,RN17B,RN17C,RN18,RN18A,RN19,RN19A,RN19B,RN20,RN20A,RN20B,
     4   BND3,RN21,RN22,RN23,RN23A,RN23B,RN25,RN25A(31),RN30A,RN30B,
     5   RN30C,RN31,BETA,RN32
      DIMENSION A1(31),A2(31)
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
        GRVT=980.
       TCA=2.43E3
       DWV=.226
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
      C76=7.66
      C358=35.86
      C172=17.26939
      C409=4098.026
      C218=21.87456
      C580=5807.695
      C610=6.1078E3
      C149=1.496286E-5
      C879=8.794142
      C141=1.4144354E7
C***   DEFINE THE COEFFICIENTS USED IN TERMINAL VELOCITY
       IF (LIN .EQ. 1) THEN
        AG=1400.
        BG=.5
       ELSE
        AG=351.2
        BG=.37
       ENDIF
CRH    AG=351.2
CRH    BG=.37
        AS=152.93
        BS=.25
CYOSH   AS=4.0858E2
CYOSH   BS=.57
CRH     AS=68.63154
CRH     BS=.11
         AW=2115.
         BW=.8
c
       if (ice913 .eq. 1 .AND. LIN .EQ. 0) then
          ag=372.3
          bs=.11
       endif
c
       BGH=.5*BG
       BSH=.5*BS
       BWH=.5*BW
       BGQ=.25*BG
       BSQ=.25*BS
       BWQ=.25*BW
C***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
      ROQR=1.
      TNW=.08
CLORD TNW=.22
        ROQS=.084
CYOSH   ROQS=.084
        TNS=.0765
CRH     TNS=.04
CYOSH   TNS=.0765
CFRED   TNS=1.
c
CRH       ROQG=.4
CYOSH     ROQG=.9
CLORD     ROQG=.3
c
CRH       TNG=.04
CFRED     TNG=.22
         IF (LIN .EQ. 1) THEN
          ROQG=.9
          TNG=.002
         ELSE
          ROQG=.4
          TNG=.1
         ENDIF
      GA3B=4.6941552
      GA4B=17.83779
      GA6B=496.6041
      GA5BH=1.827363
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
CCCCCC        LIN ET AL., 1983 OR LORD ET AL., 1984   CCCCCCCCCCCCCCCCC
      AC1=AW
      CC1=AS
      CD1=6.E-1
       CD2=4.*GRVT/(3.*CD1)
      ZRC=(CPI*ROQR*TNW)**0.25
      ZSC=(CPI*ROQS*TNS)**0.25
      ZGC=(CPI*ROQG*TNG)**0.25
      VRC=AC1*GA4B/(6.*ZRC**BW)
      VSC=CC1*GA4D/(6.*ZSC**BS)
      IF (LIN .EQ. 1) THEN
         VGC=GA4G*SQRT(CD2*ROQG/ZGC)/6.
      ELSE
         VGC=AG*GA4G/(6.*ZGC**BG)
      ENDIF
CRH   VGC=AG*GA4G/(6.*ZGC**BG)
CFRED VGC=AG*GA4G/(6.*ZGC**BG)
C     ****************************
      RN1=1.E-3
      RN2=1.E-3
      IF (LIN .EQ. 1) THEN
       BND1=7.5E-4
       BND2=1.25E-3
      ELSE
       BND1=1.E-3
       BND2=2.00E-3
c
       if (ice913 .eq. 1)  bnd2=1.5e-3
c
      ENDIF
CYOSH  BND2=3.E-3

      RN3=.25*CPI*TNS*CC1*GA3D
       ESW=1.
      RN4=.25*CPI*ESW*TNS*CC1*GA3D
CLIN       ERI=1.
       ERI=0.1
      RN5=.25*CPI*ERI*TNW*AC1*GA3B
       AMI=1./(24.*4.19E-10)
      RN6=CPI2*ERI*TNW*AC1*ROQR*GA6B*AMI
       ESR=1.
ctao
       ESR=0.5*ESR
      RN7=CPI2*ESR*TNW*TNS*ROQS
       ESR=1.
      RN8=CPI2*ESR*TNW*TNS*ROQR
       ESR=0.1*ESR
      RN9=CPI2*TNS*TNG*ROQS
      RN10=2.*CPI*TNS
       RN101=.31*GA5DH*SQRT(CC1)
       RN10A=ALS*ALS/RW
      RN11=2.*CPI*TNS/ALF
       RN11A=CW/ALF
CLIN   AMI50=4.8E-7
       AMI50=4.8E-7*(100./50.)**3
       AMI40=3.84E-9
       EIW=1.
       UI50=100.
CLIN   RI50=5.E-3
       RI50=2.*5.E-3
       CMN=1.05E-15
      RN12=CPI*EIW*UI50*RI50**2
      DO 10 K=1,31
       Y1=1.-A2(K)
      RN13(K)=A1(K)*Y1/(AMI50**Y1-AMI40**Y1)
       RN12A(K)=RN13(K)/AMI50
       RN12B(K)=A1(K)*AMI50**A2(K)
       RN25A(K)=A1(K)*CMN**A2(K)
   10 CONTINUE
       EGW=1.
      RN14=.25*CPI*EGW*TNG*GA3G*SQRT(CD2*ROQG)
CLIN    EGI=.1
        EGI=0.05
      RN15=.25*CPI*EGI*TNG*GA3G*SQRT(CD2*ROQG)
CLIN   EGI=1.
       EGI=0.5
      RN15A=.25*CPI*EGI*TNG*GA3G*SQRT(CD2*ROQG)
       EGR=1.
      RN16=CPI2*EGR*TNG*TNW*ROQR
      RN17=2.*CPI*TNG
       RN17A=.31*GA5GH*(CD2*ROQG)**.25
       RN17B=CW-CI
       RN17C=CW
       APRI=.66
       BPRI=1.E-4
      RN18=20.*CPI2*BPRI*TNW*ROQR
       RN18A=APRI
      RN19=2.*CPI*TNG/ALF
       RN19A=.31*GA5GH*(CD2*ROQG)**.25
       RN19B=CW/ALF
      RN20=2.*CPI*TNG
       RN20A=ALS*ALS/RW
       RN20B=.31*GA5GH*(CD2*ROQG)**.25
       BND3=2.E-3
      RN21=1.E3*1.569E-12/0.15
       ERW=1.
      RN22=.25*CPI*ERW*AC1*TNW*GA3B
      RN23=2.*CPI*TNW
       RN23A=.31*GA5BH*SQRT(AC1)
       RN23B=ALV*ALV/RW
       CN0=1.E-8
      RN25=CN0
       RN30A=ALV*ALS*AMW/(TCA*ARS)
       RN30B=ALV/TCA
       RN30C=ARS/(DWV*AMW)
      RN31=1.E-17
       BETA=-.6
      RN32=4.*51.545E-4
c     ****************************
      do 20 k=1,nz
       srr(k)=1./sqrt(rho(k))
   20  qrr(k)=sqrt(srr(k))
      RETURN
      END
