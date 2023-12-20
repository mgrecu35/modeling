c-----------------------------------------------------------------------
      subroutine obs_toga(sst,ssp,press,temp,vap,uu,vv,ww,ubo_6h,vbo_6h,
     1                    wbo_6h,q1o_6h,q2o_6h)
      parameter(NPIN=41,ITT=244)
      parameter(ISKIPNOV=1)
      parameter(ISKIPDEC=0)
      parameter(ISKIPJAN=0)
      parameter(ISKIPDAY=18)
      parameter(ISTARTHOUR=0)
      dimension ttt(npin,itt),qqq(npin,itt),ppp(npin,itt),sst(itt),
     1          ssp(itt)
      dimension press(npin),temp(npin),vap(npin),uu(npin),vv(npin),
     1          ww(npin),ubo_6h(npin,itt),vbo_6h(npin,itt),
     2          wbo_6h(npin,itt),q1o_6h(npin,itt),q2o_6h(npin,itt)
      common/upper/ t_adjust(npin,itt),q_adjust(npin,itt),press_t(npin),
     1      temp_t(npin),vap_t(npin)
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      open(60,file='/usr/raid1/djohnson/toga.ifa.480days/fields.ifa'
     1       ,status='old')
      open(61,file='/usr/raid1/djohnson/toga.ifa.480days/advect.ifa'
     1       ,status='old')
      open(62,file='/usr/raid1/djohnson/toga.ifa.480days/vert.force'
     1       ,status='old')
      open(63,file='/usr/raid1/djohnson/toga.ifa.480days/misc.ifa'
     1       ,status='old')

C  Read into arrays, the times at beginning of simulation to be written over

      iskip = (iskipnov*30 + iskipdec*31 + iskipjan*31 + iskipday) * 4
     1      + istarthour/6


      do i=1,iskip
         read(60,*)
         read(61,*)
         read(62,*)
         read(63,*)
         do k=1,npin
            read(60,*)
            read(61,*)
            read(62,*)
         enddo
      enddo

C  Now read in the data we want to use in the simulation


      do 999 i=1,itt

         read(60,*)
         read(61,*)
         read(62,*)
         read(63,*) dum,dum,dum,dum,dum,sst(i)

         do k=1,npin

            read(60,*) ppp(k,i),dum,ubo_6h(k,i),vbo_6h(k,i),wbo_6h(k,i),
     1                 ttt(k,i),dum,qqq(k,i)
            read(61,*) dum,dum,dum,dum,dum,ht,dum,hq,vq
            read(62,*) dum,dum,dum,vt
            q1o_6h(k,i) = (ht + vt) * 86400.
            q2o_6h(k,i) = (hq + vq) * 86400. * 1000.
         enddo

         PRINT*,'NON-MODIFIED LARGE SCALE CONDITIONS AT',I*6-6,'HOURS'
         WRITE(6,240)
         DO K1=1,NPIN-1
            K=NPIN+1-K1
            WRITE(6,200) K,PPP(K,I),TTT(K,I),QQQ(K,I),UBO_6H(K,I),
     2                   VBO_6H(K,I),WBO_6H(K,I),Q1O_6H(K,I),Q2O_6H(K,I)
         ENDDO

        ssp(i) = ppp(1,i)
        PSFC=SSP(1)*1000.

        TTT(1,I)=((TTT(2,I)+273.16)*(PPP(1,I)/PPP(2,I))**(RA/CP))-273.16
        IF (TTT(NPIN-2,I)-TTT(NPIN-3,I) .GT. 10.)
     1                                   TTT(NPIN-2,I)=TTT(NPIN-3,I)+10.
        TTT(NPIN-1,I)=TTT(NPIN-2,I)+10.
        TTT(NPIN,I)=TTT(NPIN-1,I)+10.

        QQQ(1,I)=QQQ(2,I)
        QQQ(NPIN-1,I)=0.
        QQQ(NPIN,I)=0.

        DO K=35,NPIN-3
          IF (UBO_6H(K,I) .LE. -20.0) UBO_6H(K,I)=-20.0
        ENDDO

        DO K=2,NPIN-2
          IF (PPP(K,I) .Le. 100.) Q1O_6H(K,I)=0.
          IF (PPP(K,I) .Le. 100.) Q2O_6H(K,I)=0.
          IF (PPP(K,I) .Le. 100.) QQQ(K,I)=0.
          IF (PPP(K,I) .Le. 100.) WBO_6H(K,I)=0.0
          if (i .ge. 5) then
            IF (PPP(K,I) .LE. 125.) THEN
c             WBO_6H(K,I)=0.75*((REAL(NPIN)-REAL(K))/9.)**2*WBO_6H(K,I)
             Q1O_6H(K,I)=0.5*Q1O_6H(K,I)
             Q2O_6H(K,I)=0.5*Q2O_6H(K,I)
             if (Q1O_6H(K,I) .ge. 1.5) Q1O_6H(K,I)=1.5
             if (Q1O_6H(K,I) .le. -1.5) Q1O_6H(K,I)=-1.5
            ENDIF
          else
            IF (PPP(K,I) .LE. 150.) THEN
             WBO_6H(K,I)=0.
             Q1O_6H(K,I)=0.
             Q2O_6H(K,I)=0.
            ENDIF
          endif
        ENDDO

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        dpp=(ppp(1,i)-ppp(2,i))/(ppp(2,i)-ppp(3,i))
        ubo_6h(1,i)=ubo_6h(2,i)+(ubo_6h(2,i)-ubo_6h(3,i))*dpp
        if ( abs(ubo_6h(1,i)) .ge. abs(ubo_6h(2,i)) )
     1      ubo_6h(1,i)=ubo_6h(2,i)
        ubo_6h(npin-3,i)=ubo_6h(npin-4,i)
        ubo_6h(npin-2,i)=ubo_6h(npin-3,i)
        ubo_6h(npin-1,i)=ubo_6h(npin-2,i)
        ubo_6h(npin,i)=ubo_6h(npin-1,i)

        vbo_6h(1,i)=vbo_6h(2,i)+(vbo_6h(2,i)-vbo_6h(3,i))*dpp
        if ( abs(vbo_6h(1,i)) .ge. abs(vbo_6h(2,i)) )
     1     vbo_6h(1,i)=vbo_6h(2,i)
        vbo_6h(npin-3,i)=vbo_6h(npin-4,i)
        vbo_6h(npin-2,i)=vbo_6h(npin-3,i)
        vbo_6h(npin-1,i)=vbo_6h(npin-2,i)
        vbo_6h(npin,i)=vbo_6h(npin-1,i)

        wbo_6h(1,i)=0.
        wbo_6h(npin-1,i)=0.
        wbo_6h(npin,i)=0.

        q1o_6h(1,i)=q1o_6h(2,i)+(q1o_6h(2,i)-q1o_6h(3,i))*dpp
        if ( abs(q1o_6h(1,i)) .ge. abs(q1o_6h(2,i)) )
     1     q1o_6h(1,i)=q1o_6h(2,i)
        q1o_6h(npin-1,i)=0.
        q1o_6h(npin,i)=0.

        q2o_6h(1,i)=q2o_6h(2,i)+(q2o_6h(2,i)-q2o_6h(3,i))*dpp
        if ( abs(q2o_6h(1,i)) .ge. abs(q2o_6h(2,i)) )
     1     q2o_6h(1,i)=q2o_6h(2,i)
        q2o_6h(npin-1,i)=0.
        q2o_6h(npin,i)=0.

999   continue

      DO K=1,NPIN
        PRESS_T(K)=0.
        TEMP_T(K)=0.
        VAP_T(K)=0.
      ENDDO
      SSP_T=0.
      SST_T=0.

      DO I=1,ITT
        IP=I+1
        IF(I.EQ.ITT) IP=ITT
        SSP_T=SSP_T+SSP(I)
        SST_T=SST_T+SST(I)
        DO K=1,NPIN
          PRESS_T(K)=PRESS_T(K)+PPP(K,I)
          TEMP_T(K)=TEMP_T(K)+TTT(K,I)
          VAP_T(K)=VAP_T(K)+QQQ(K,I)
          T_ADJUST(K,I)=TTT(K,IP)-TTT(K,I)
          Q_ADJUST(K,I)=QQQ(K,IP)-QQQ(K,I)
        ENDDO
      ENDDO

      DO K=1,NPIN
        PRESS_T(K)=PRESS_T(K)/ITT
        TEMP_T(K)=TEMP_T(K)/ITT
        VAP_T(K)=VAP_T(K)/ITT
      ENDDO
      SSP_T=SSP_T/ITT
      SST_T=SST_T/ITT

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
  240 format(//,'level          p           t           q           u
     1       v           w         q1        q2
     2 ')
  241 format(//,'level          p           t           q           u
     1       v           w')
      end
