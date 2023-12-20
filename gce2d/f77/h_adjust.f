       
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
