Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bulk_flux(hUm,hTm,hUs,hTs,ws,sst,atb,qq,ws_h,Ta_h,qq_h,
     1                     rs,rl,rainx,pp,zix,Jcoolx,Jwarmx,HF,EF,RF,
     2                     TAU,Ustar,Tstar,Qstar,CD,CH,CE,RRx,RTx,RQx,
     3                     ZLx,ZOx,zotx,zoqx,dt_wrmx,dterx,T0)

      real hUm,hTm,hUs,hTs,ws_h,Ta_h,qq_h,ws,sst,atb,qq,pp,zix,rainx
      real HF,EF,TAU,Ustar,Qstar,Tstar,rl,rs,RF,T0,CD,CE,CH,RRx,RTx,RQx
      real Zlx,ZOx,Jcoolx,Jwarmx,zotx,zoqx,dt_wrmx,dterx
      real jtime,qcol_ac,tau_ac
      integer jamset
      common /old/jtime,qcol_ac,tau_ac,tau_old,rf_old,hf_old,ef_old
     1            ,jamset
      COMMON/PIN/U,T,Q,TS,QS,rns,rnl,ZU,ZT,ZQ,zi,P,ID
      COMMON/POUT/USR,TSR,QSR,ZO,zot,zoq,ZL,RR,RT,RQ,RI,dter,dqer,tkt
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
      common/wgust/DU_Wg,Wg     !MC

      common/kier/IER          ! cccshie
      integer IER          ! cccshie

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
c -------------------- correct SST with PWP model ------------
      if(jwarm.eq.0) go to 15                   ! by-pass warm layer
Cwang if(jwarm.eq.2.) then                      ! first line of data
Cwang    jump=1                                 
Cwang    go to 16                               ! set jtime and pass thru' ASL
Cwang end if
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
Cwang&     tk_pwp=min(19.,ctd1*tau_ac/sqrt(qcol_ac+qjoule))
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
c        CALL ASL(Jcool,IER)
         CALL ASL(Jcool) ! shie 5/7/01, move "ier" to common block
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
