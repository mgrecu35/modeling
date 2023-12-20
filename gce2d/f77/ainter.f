C-----------------------------------------------------------------------
      subroutine ainter(sst,ssp,tb,fd,qb,fe,ub,vb,wb,y1,y2,toptq,topuv)
      PARAMETER (N=41,ITT=244,NZ=43,KLES=NZ-1,KMAX1=KLES+1)
      PARAMETER (KNN=18*N,NN1=N+1,NN2=2*N+1,N4=N*4)
      common/msui/ isounding,isui,ifine,idelz
      dimension y1(nz),y2(nz),itab(3),iop(2)
      real fd(nz),fe(nz),tb(nz),qb(nz),ub(nz),vb(nz),wb(nz)
      real uo_6h(n,itt),vo_6h(n,itt),wo_6h(n,itt),q1o_6h(n,itt)
      real q2o_6h(n,itt)
      real P(N),T(N),Q(N),U(N),V(N),WW1(N),H1(N),H2(N)
      real TAB(3),rqo,AT(162),AQ(162),H(N),WK(N4),W(N)
      real Z(nz),Z1(nz),ZZ(168)
      real QQ1(KNN)
      common/bb6/ tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     $   qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     $   ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt)
     $   ,q2t(nz),vb_6h(nz,itt),vbt(nz)
      COMMON/UPPER/ T_ADJUST(N,ITT),Q_ADJUST(N,ITT),P_T(N),
     1      T_T(N),V_T(N)
      common/upper1/ t_adjust0(nz,itt),q_adjust0(nz,itt),press_t0(nz),
     1     temp_t0(nz),vap_t0(nz)
      common /zzzobs/ zzz(n,itt)
      real sst(itt),ssp(itt),dum(n)
      EQUIVALENCE (U(1),QQ1(1)),(V(1),QQ1(NN1)),(WW1(1),QQ1(NN2))
    
      character*4 IVAA(9)
      DATA IVAA/'P  ',' H ',' T ',' Q ',' U ',' V ','WW1 ','RHW '
     1          ,'RHI '/
      save
      call obs(sst,ssp,p,t,q,u,v,ww1,uo_6h,vo_6h,wo_6h,
     1              q1o_6h,q2o_6h)
c      call obs_toga(sst,ssp,p,t,q,u,v,ww1,uo_6h,vo_6h,wo_6h,
c     1              q1o_6h,q2o_6h)
      KMAX=KLES+KLES+1
      DO K=1,KLES
        Z(K)=y2(k+1)*.01
        Z1(K)=y1(k+1)*.01
      ENDDO
      Z(KMAX1)=Z(KLES)+(Z(KLES)-Z(KLES-1))
      Z1(KMAX1)=Z1(KLES)+(Z1(KLES)-Z1(KLES-1))
      write(6,2031)
      DO 11 K=1,KMAX1
   11   write(6,100) K,Z(K),Z1(K)
      DO 12 K=1,KLES
        I=K*2
        I1=K*2-1
        ZZ(I1)=Z(K)
   12   ZZ(I)=Z1(K)
      ZZ(KMAX)=Z(KMAX1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      h(1)=0.
      term1=2.87e+6/980.*.5
      do k=2,n
      km=k-1
      h(k)=h(km)-term1*((t_t(k)+273.16)*(1.+.61*v_t(k)*1.e-3)
     1   +(t_t(km)+273.16)*(1.+.61*v_t(km)*1.e-3))*log(p_t(k)/p_t(km))
     2     *.01
       add=3.799e3/p_t(km)*exp(17.26939-4098.026/(t_t(km)+273.16-35.86))
      h1(km)=v_t(km)/add
       add1=3.799e3/p_t(km)*exp(21.87456-5807.695/(t_t(km)+273.16-7.66))
      h2(km)=v_t(km)/add1
      enddo
      write(6,*)'      P           H       T      Q       H1      H2'
      do k=1,n
        write(6,1031)k,p(k),h(k),t(k),q(k),h1(k),h2(k)
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      H(1)=0.
      TERM1=2.87E+6/980.*.5
cccshie 4/9/01 using initial time obs height zzz(k,1)
      DO 20 K=2,N
      KM=K-1
c     H(K)=H(KM)-TERM1*((T(K)+273.16)*(1.+.61*Q(K)*1.E-3)
c    1     +(T(KM)+273.16)*(1.+.61*Q(KM)*1.E-3))*LOG(P(K)/P(KM))*.01
c      print*,'h,z',h(k),zzz(k,1)
      h(k)=zzz(k,1)  ! 4/9/01
       ADD=3.799E3/P(KM)*EXP(17.26939-4098.026/(T(KM)+273.16-35.86))
      H1(KM)=Q(KM)/ADD
       ADD1=3.799E3/P(KM)*EXP(21.87456-5807.695/(T(KM)+273.16-7.66))
      H2(KM)=Q(KM)/ADD1
   20 CONTINUE
      zpresst=100.
      zpressu=100.
      do k=1,n
        if(p(k).ge.zpresst) toptq=h(k)
        if(p(k).ge.zpressu) topuv=h(k)
      enddo
      print*,'top at or just below 100mb =',toptq
      print*,'top at or just below 100mb =',topuv
      write(6,102) (IVAA(I),I=1,9)
      DO 30 K=1,N
   30 write(6,1031)K,P(K),H(K),T(K),Q(K),U(K),V(K),WW1(K),H1(K),H2(K)
      do k=1,n
        if(h(k).le.z(nz-1)) ntop=k
      enddo
      rqo=0.
      do k=1,n-1
        rqo=rqo+.5*(q(k)+q(k+1))*.001*(p(k)-p(k+1))*1000./980.
      enddo
       print*
       print*,'OBSERVED QV INTEGRATION=',rqo
      IOP(1)=4
      IOP(2)=4
      INT=1
      ITAB(1)=1
      ITAB(2)=0
      ITAB(3)=0
      CALL COEFF (N,H,T,W,IOP,INT,WK)
      DO K=1,KMAX
        Y=ZZ(K)
        CALL TERP1 (N,H,T,W,Y,INT,TAB,ITAB)
        AT(K)=TAB(1)
      ENDDO
      CALL COEFF (N,H,Q,W,IOP,INT,WK)
      write(6,1021) (IVAA(I),I=2,4)
      DO 250 K=1,KMAX
        Y=ZZ(K)
       CALL TERP1 (N,H,Q,W,Y,INT,TAB,ITAB)
       AQ(K)=TAB(1)
        if(y .ge. 16000.) aq(k)=0.
       write(6,106) K,Y,AT(K),AQ(K)
  250 CONTINUE
      DO KK=1,KLES
         kkp=kk+1
        KK1=(KK-1)*2+1
        KK2=KK1+1
       FD(KKp)=AT(KK1)
       TB(KKp)=AT(KK2)
       FE(KKp)=AQ(KK1)
       QB(KKp)=AQ(KK2)
       IF (QB(KKp) .LE. 0.0) QB(KKp)=0.
       IF (FE(KKp) .LE. 0.0) FE(KKp)=0.
      ENDDO

      CALL COEFF (N,H,U,W,IOP,INT,WK)
      do K=1,KLES
         Y=Z1(K)
         CALL TERP1 (N,H,U,W,Y,INT,TAB,ITAB)
         UB(K+1)=TAB(1)*100.
      enddo

      CALL COEFF (N,H,V,W,IOP,INT,WK)
      do K=1,KLES
         Y=Z1(K)
         CALL TERP1 (N,H,V,W,Y,INT,TAB,ITAB)
         VB(K+1)=TAB(1)*100.
      enddo

      CALL COEFF (N,H,WW1,W,IOP,INT,WK)
      do K=1,KLES
         Y=Z(K)
         CALL TERP1 (N,H,WW1,W,Y,INT,TAB,ITAB)
         WB(K+1)=TAB(1)
      enddo



      WRITE(6,1022) IVAA(2),(IVAA(I),I=5,6),IVAA(2),IVAA(7)
      DO K=1,KMAX1
         WRITE(6,103) K,Z1(K),UB(K),VB(K),Z(K),WB(K)
      ENDDO


      DO 700 KK=1,ITT
       do k=1,n
         dum(k)=uo_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 73 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        UB_6H(K+1,kk)=TAB(1)*100.
   73 CONTINUE
       do k=1,n
         dum(k)=vo_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 83 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        VB_6H(K+1,kk)=TAB(1)*100.
   83 CONTINUE
       do k=1,n
         dum(k)=wo_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 86 K=1,KLES
        Y=Z(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        WB_6H(K+1,kk)=TAB(1)
   86 CONTINUE
       do k=1,n
         dum(k)=q1o_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 88 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        Q1_6H(K+1,kk)=TAB(1)
   88 CONTINUE
       do k=1,n
         dum(k)=q2o_6h(k,kk)
       enddo
      CALL COEFF(N,H,DUM,W,IOP,INT,WK)
      DO 89 K=1,KLES
        Y=Z1(K)
        CALL TERP1(N,H,DUM,W,Y,INT,TAB,ITAB)
        Q2_6H(K+1,kk)=TAB(1)*.001
   89 CONTINUE

       do k=1,n
         dum(k)=t_adjust(k,kk)
       enddo
      call coeff(n,h,dum,w,iop,int,wk)
      do 173 k=1,kles
        y=z1(k)
        call terp1(n,h,dum,w,y,int,tab,itab)
        t_adjust0(k+1,kk)=tab(1)
  173 continue

       do k=1,n
         dum(k)=q_adjust(k,kk)
       enddo
      call coeff(n,h,dum,w,iop,int,wk)
      do 183 k=1,kles
        y=z1(k)
        call terp1(n,h,dum,w,y,int,tab,itab)
        q_adjust0(k+1,kk)=tab(1)
  183 continue

  700 CONTINUE

      RETURN
  100 FORMAT(6X,I5,2F10.2,2x,F10.8,2x,F10.8)
  102 FORMAT(//,5X,6A8,3A7)
 1021 FORMAT(//,6X,A8,A11,A10)
 1022 FORMAT(//,6X,a8,2a10,a8,A10)
  103 FORMAT(1X,I2,12F10.2)
 1031 FORMAT(I3,F9.2,F10.2,F8.2,F7.2,2F8.2,F7.2,2F6.2)
  106 FORMAT(2X,I4,8F10.3)
2031   format('                          in ainter',/,'          k     
     *z(k)     z1(k)')
      END
