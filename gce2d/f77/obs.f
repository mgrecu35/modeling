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
c     parameter(ISTARTDAY=16) !  4/3/01 shie, skip 15 days for 1st 10days
c     parameter(ISTARTDAY=18) ! 8/1/01 ,skip 17 days for 1st 2days 5/18-5/20
c     parameter(ISTARTDAY=17) ! 6/7/01, skip 16 days for a 5-day run 5/17-5/22
c     parameter(ISTARTDAY=20) ! 6/29/01, 2nd 4-5 day run 5/17-5/22
c     parameter(ISTARTDAY=33) ! 6/7/01, skip 32 days for a 9-day run 6/2-6/11
c     parameter(ISTARTDAY=14) ! 5/3/01, may14.12z-may24.12z, NX 258
c     parameter(ISTARTDAY=20) ! 4/18/01,skip 19 days for 1st 10days
c     parameter(ISTARTDAY=18) !  4/9/01 shie, skip 17 days for 1st 5days
c     parameter(ISTARTDAY=21) !  4/18/01 shie, skip 20 days for 1st 10days
c     parameter(ISTARTDAY=21) !  4/8/01 shie, skip 20 days for 1st 5days
c     parameter(ISTARTDAY=16) ! 3/16/01 shie, skip 15 days for 11-20 days 
c     parameter(ISTARTHOUR=0) ! 4/20, 4/18, 4/8, 6/7
c     parameter(ISTARTHOUR=6) ! 4/20, 4/18, 4/8, 6/7
c     parameter(ISTARTHOUR=12) ! 4/4 ! 5/3/01,may14.12z-may24.12z, NX 258
c     parameter(ISTARTHOUR=18) ! 4/5
c     parameter(ISTARTHOUR=8)  ! 4/6 move forward 8 hrs
c     parameter(ISTARTHOUR=8)  ! 4/9, 4/18 move forward 8 hrs
c     parameter(ISTARTHOUR=14)  ! 4/18 move forward 14 hrs

c     parameter(ISTARTDAY=57)  ! 8/21/01, 062612-062900 (31+26=57)
c     parameter(ISTARTDAY=44)  ! 8/21/01, 061312-061412 (31+13=44)
c     parameter(ISTARTDAY=14)  ! 8/21/01, 051412-051512 ( 0+14=14)
c     parameter(ISTARTDAY=27)  ! 8/21/01, 052712-052812 ( 0+27=27)
c     parameter(ISTARTHOUR=12) ! 8/21/01, 062612 (12z)
      parameter(ISTARTDAY=18)  ! 8/29/01, 051800-052600 ( 0+18=18)
c     parameter(ISTARTDAY=33)  ! 8/29/01, 060200-061100 (31+ 2=33)
      parameter(ISTARTHOUR=0)  ! 8/29/01, 051800-052600 ( 0+18=18)

cccshie average over the troble data at 052018(may 20, 18z):  iavge=1
cc     iavge=1  ! average over the trouble data at 052018 (old data only)
       iavge=0  ! use original trouble data at 052018, or use new data 6/5/01

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

      DO I=1,ITT
        IP=I+1
        IF(I.EQ.ITT) IP=ITT
        DO K=1,NPIN
          PRESS_T(K)=PRESS_T(K)+PPP(K,I)
          TEMP_T(K)=TEMP_T(K)+TTT(K,I)
          VAP_T(K)=VAP_T(K)+QQQ(K,I)
        ENDDO
      ENDDO

      DO K=1,NPIN
        PRESS_T(K)=PRESS_T(K)/ITT
        TEMP_T(K)=TEMP_T(K)/ITT
        VAP_T(K)=VAP_T(K)/ITT
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
