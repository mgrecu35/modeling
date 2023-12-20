c ---------------------------------------------------------------------
c     SUBROUTINE ASL(Jcool,IER)
      SUBROUTINE ASL(Jcool)   ! cccshie
c
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
      common/wgust/DU_Wg,Wg      ! MC

      common/kier/IER          ! cccshie
      integer IER          ! cccshie

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
       IF(RI.GT.0.25) IER=-1
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
