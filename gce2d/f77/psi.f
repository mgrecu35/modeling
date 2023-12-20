c
c
      FUNCTION PSI(ID,ZL)
C
C       TO EVALUATE THE STABILITY FUNCTION PSI FOR WIND SPEED (IFLAG=1)
C       OR FOR TEMPERATURE AND HUMIDITY PROFILES FROM STABILITY 
C       PARAMETER ZL. SEE LIU ET AL (1979).
c       Modified to include convective form following Fairall (Unpublished)
C
      IF(ZL)10,20,30                                                 
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
