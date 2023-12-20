Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ZETA(T,Q,USR,TSR,QSR,Z,ZL)
C
C       TO EVALUATE OBUKHOVS STABILITY PARAMETER Z/L FROM AVERAGE
C       TEMP T IN DEG C, AVERAGE HUMIDITY Q IN GM/GM, HEIGHT IN M,
C       AND FRICTIONAL VEL,TEMP.,HUM. IN MKS UNITS
C       SEE LIU ET AL. (1979)
C
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
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
