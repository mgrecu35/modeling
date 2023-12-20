Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c      subroutine SFflux (mmmm,ATIME,UBOT,VBOT,HT,SSTIN,TBOT,
      subroutine SFflux (mmmm,UBOT,VBOT,HT,SSTIN,TBOT,
     1                   QBOT,RSHRT,RLONG,RAINR,PSFC)
      PARAMETER (NX=514,NZ=43)
      real    p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),pairsfc(nx)
     1,        thairsf(nx)
      common/dinrad/ p00,dz0,tairsfc,qairsfc,thairsf,pairsfc

      COMMON /SFLUXS/ SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx
      DIMENSION UBOT(NX),VBOT(NX),TBOT(NX),QBOT(NX),RAINR(NX)  
      DIMENSION RSHRT(NX),RLONG(NX),SSTIN(NX)
      DIMENSION USTARs(NX),TSTARs(NX),QSTARs(NX)
      common/kier/IER          ! cccshie
      integer IER          ! cccshie

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

      real hUm,hTm,hUs,hTs,ws_h,qq_h,ta_h
c      real loc
      real ws,sst,atb,qq,pp,zi,rain
      real QH,QE,TAU,Ustar,Qstar,Tstar
      real rl,rs,RF,T0
c      real ts_depth
      real CD,CE,CH,RR,RT,RQ,Zl,ZO
      real Jcool,Jwarm,zot,zoq,dt_wrm,dter
c      real time
      real jtime,qcol_ac,tau_ac
      integer jamset
      common /old/jtime,qcol_ac,tau_ac,tau_old,rf_old,hf_old,ef_old
     1,           jamset
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
c      time=ATIME    ! fraction of a day (<=1)
      zi=600.       ! m
      pp=PSFC       ! mb

      aveqh=0.
      aveqe=0.
      aveustar=0.

      fnum=0.0  ! cccshie
      do i=2,mmmm
       ws= SQRT(UBOT(I)**2+VBOT(I)**2)      ! m/s
       sst=SSTIN(i)                         ! degrees C
       atb=TBOT(I)                          ! degrees C
       qq=QBOT(I)                           ! g/kg
       rs=RSHRT(I)                          ! W/m^2
       rl=RLONG(I)                          ! W/m^2
       rain=RAINR(I)                        ! mm/hr
c1      read(10,*,end=999) time,ws,sst,atb,qq,rs,rl,rain,pp,zi
       call bulk_flux(hUm,hTm,hUs,hTs,ws,sst,atb,qq,WS_H,TA_H,QQ_H,rs,
     1                rl,rain,pp,zi,jcool,jwarm,QH,QE,RF,TAU,USTAR,
     2                TSTAR,QSTAR,CD,CH,CE,RR,RT,RQ,ZL,ZO,ZOT,ZOQ,
     3                DT_WRM,DTER,T0)

cccshie
       if(ier.ge. 0) then
        fnum=fnum+1.0
       aveqh=aveqh+qh
       aveqe=aveqe+qe
       aveustar=aveustar+ustar
       endif

       ustarS(I)=ustar
       tstarS(I)=tstar
       qstarS(I)=qstar
      enddo

       if(fnum.ne.0.) then
       aveqh=aveqh/fnum
       aveqe=aveqe/fnum
       aveustar=aveustar/fnum
       else
       aveqh=-999.
       aveqe=-999.
       aveustar=-999.
       endif

	print *, 'fnum, sensible heat flux, latent heat flux, ustar'
       print *, fnum,aveqh,aveqe,aveustar

c      aveqh=aveqh/float(nx-2)
c      aveqe=aveqe/float(nx-2)
c      aveustar=aveustar/float(nx-2)
c      print *, 'sensible heat flux, latent heat flux, and ustar'
c      print *, aveqh, aveqe,aveustar

C
      do i=1,nx
       suw(i)=0.
       svw(i)=0.
       swt(i)=0.
       swq(i)=0.
      enddo
C
c       ws= SQRT(UBOT(2)**2+VBOT(2)**2)            
c       suw(2)=-(2.*ustars(2)-ustars(3))**2*ubot(2)/ws 
c       svw(2)=-(2.*ustars(2)-ustars(3))**2*vbot(2)/ws
c       SWT(2)=-USTARS(2)*TSTARS(2)
c       SWQ(I)=-USTARS(2)*QSTARS(2)


      do i=2,mmmm
        ws= SQRT(UBOT(I)**2+VBOT(I)**2)
cccshie
	if(ustars(i).ne.-999. .and. TSTARS(I).ne.-999. 
     1   .and. QSTARS(I).ne.-999.) then
       SUW(I)=-ustars(i)**2*UBOT(I)/WS
       SVW(I)=-ustars(i)**2*VBOT(I)/WS
       SWT(I)=-USTARS(I)*TSTARS(I)
       SWQ(I)=-USTARS(I)*QSTARS(I)
        endif
      enddo
      
c      write(29) suw
c      write(29) svw
c      write(29) swt
c      write(29) swq

      return
      end
