cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE PBLIN (IFLAG,isfcave)
C     ****   SET DOMAIN AND INITIAL CONDITION
      PARAMETER (NX=514,NZ=43,iles=nx-1)

      real dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      real    sec,aminut,rdt
      common/rbstart/ sec,aminut,rdt
      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1
      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),pairsfc(nx)
     1,        thairsf(nx)
      common/dinrad/ p00,dz0,tairsfc,qairsfc,thairsf,pairsfc

      real ri180(nx),riold180(nx)
      integer iricont
      common/sfcri/ ri180,riold180,iricont

      common/dinrad1/ dz1half
      COMMON/BPBL/ UHT(NZ),WHT(NZ),TGBAT0

      COMMON/B2U/ U(NX,NZ)
      COMMON/B2V/ V(NX,NZ)
      COMMON/B1T/ DPT(NX,NZ)
      COMMON/B1Q/ DQV(NX,NZ)
      COMMON/B1C/ QCL(NX,NZ)
      COMMON/B1I/ QCI(NX,NZ)
      real    SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      COMMON/SFLUXS/ SUW,SVW,SWT,SWQ

      DIMENSION UBOT(NX),VBOT(NX),TBOT(NX),QBOT(NX),RAINR(NX)
      DIMENSION RSHRT(NX),RLONG(NX),SSTIN(NX)

      IF (IFLAG.EQ.0) RETURN

         arii=0.
         do i=2,iles
          arii=arii+ri180(i)
         enddo
          arii=arii*ril2
         mmmm=iles
          if (isfcave .eq. 1) mmmm=2
       DO  I=2,mmmm
        if (isfcave .eq. 0) then
          ip=i+1
           if (i .eq. iles) ip=2
          UBOT(I)=.005*(U(I,2)+u(ip,2))       !to m/s
          VBOT(I)=.005*(V(I,2)+v(ip,2))
          TBOT(I)=(TA1(2)+DPT(I,2))*pi(2)-273.16      ! potential T
          QBOT(I)=1000.*(QA1(2)+DQV(I,2))     ! vapor mixing ratio (g/kg)
          RSHRT(I)=140.
          RLONG(I)=400.                       !for now
          RAINR(I)=ri180(I)
          SSTIN(I)=tairsfc(i)-273.16
        else
          ubot(i)=.01*ub1(2)         !to m/s
          vbot(i)=.01*vb1(2)
          tbot(i)=ta1(2)*pi(2)-273.16      ! potential t
          qbot(i)=1000.*qa1(2)       ! vapor mixing ratio (g/kg)
          rshrt(i)=140.
          rlong(i)=400.              !for now
          rainr(i)=arii
          sstin(i)=tairsfc(i)-273.16
        endif
       ENDDO    
       HT=dz1half
c       PSFCMB=.001*PSFC        ! pressure at surfase (from microbar to mb)
       PSFCMB=.001*pairsfc(1)   ! pressure at surfase (from microbar to mb)
       ATIME=.2                 ! for now

c       CALL SFflux (mmmm,ATIME,UBOT,VBOT,HT,SSTIN,TBOT,
       CALL SFflux (mmmm,UBOT,VBOT,HT,SSTIN,TBOT,
     1              QBOT,RSHRT,RLONG,RAINR,PSFCMB)
       if (isfcave .eq. 1) then
         do i=3,iles
          suw(i)=suw(2)
          svw(i)=svw(2)
          swt(i)=swt(2)
          swq(i)=swq(2)
         enddo
       endif

      RETURN 
      END
