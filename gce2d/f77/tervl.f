Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TERVL (IRSG,RHO,FV)
C     COMPUTE THE TERMINAL VELOCITY OF QR, QS AND QG
      PARAMETER (NX=514,NZ=43)
      PARAMETER (NX10=10*NX)
      common/ilin/ lin

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC

      real    ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq

      real    qrn(nx,nz),qcg(nx,nz),qcs(nx,nz)
      common/b1r/ qrn
      common/b1g/ qcg
      common/b1s/ qcs

      real    ww1(nx,nz)
      COMMON/BTV/ WW1

      real    y1(nx),y2(nx),y3(nx),vr(nx),vs(nx),vg(nx),y4(nx10)
      common/ba/ y1,y2,y3,vr,vs,vg,y4

      integer irsg
      real    fv(nz),rho(nz)

C     ******************************************************************
      IF(IRSG.NE.0) GO TO 1
C     ***      IRSG=0    FOR RAIN WATER (QR)  **********
      DO 100 K=2,KLES
         VRCF=VRC*FV(K)
       DO 10 I=2,ILES
        WW1(I,K)=0.0
          A1=.5*RHO(K)*(QRN(I,K)+QRN(I,K-1))
        IF (A1 .GT. 1.E-15) THEN
         WW1(I,K)=VRCF*A1**BWQ
        ENDIF
   10  CONTINUE
  100 CONTINUE
      RETURN
CC    ***      IRSG=1    FOR SNOW (QS)        **********
    1 IF(IRSG.NE.1) GO TO 2
      DO 200 K=2,KLES
          VSCF=VSC*FV(K)
       DO 20 I=2,ILES
        WW1(I,K)=0.0
          A1=.5*RHO(K)*(QCS(I,K)+QCS(I,K-1))
         IF (A1 .GT. 1.E-15) THEN
          WW1(I,K)=VSCF*A1**BSQ
         ENDIF
   20  CONTINUE
  200 CONTINUE
      RETURN
CC    ***      IRSG=2    FOR GRAUPEL (QG)      **********
    2 DO 300 K=2,KLES
          IF (LIN .EQ. 1) THEN
            VGCR=VGC/SQRT(RHO(K))
          ELSE
            VGCR=VGC*FV(K)
          ENDIF
       DO 30 I=2,ILES
        WW1(I,K)=0.0
          A1=.5*RHO(K)*(QCG(I,K)+QCG(I,K-1))
         IF (A1 .GT. 1.E-15) THEN
          WW1(I,K)=VGCR*A1**BGQ
         ENDIF
   30  CONTINUE
  300 CONTINUE
      RETURN
      END
