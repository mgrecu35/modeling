c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sat (n,id,fv)
c     ****   compute ice phase microphysics and saturation processes
      parameter (NX=514,NZ=43,NT=2880,ITT=244)
      parameter (nt5=5*nt,nb=nx*nz-31*nx,nb1=nx*nz-7*nx)
      common/timestat/ ndt_stat,itime_ave,mmave
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      common/iice/ new_ice_sat
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/rterv/ zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      common/rsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,rn2,bnd2,rn3,rn4,
     $  rn5,rn50,rn51,rn52,rn53,rn6,rn60,rn61,rn62,rn63,rn7,rn8,rn9,
     $  rn10,rn101,rn102,rn10a,rn10b,rn10c,rn11,rn12,rn12a(31),
     $  rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn171,rn172,rn17a,rn17b,
     $  rn17c,rn18,rn18a,rn19,rn191,rn192,rn19a,rn20,rn20a,rn20b,rn30,
     $  rn30a,rn21,bnd21,rn22,rn23,rn231,rn232,rn25,rn25a(31),rn31,beta,
     $  rn32,rn33,rn331,rn332,rn34,rn35
      common/bsize/ it1,it2,kt1,kt2
      common/b1t/ dpt(nx,nz)
      common/b1q/ dqv(nx,nz)
      common/b1c/ qcl(nx,nz)
      common/b1r/ qrn(nx,nz)
      common/b2t/ dpt1(nx,nz)
      common/b2q/ dqv1(nx,nz)
      common/b2c/ qcl1(nx,nz)
      common/b2r/ qrn1(nx,nz)
      common/b2w/ ww1(nx,nz)
      common/slwave/ rsw(nx,nz),rlw(nx,nz)
      common/bsat/ pwacs(nx),wgacr(nx),dep(nx),pidep(nx),pint(nx),
     $  dd(nx),dd1(nx),qvs(nx),dm(nx),rsub1(nx),col(nx),cnd(nx),rq(nx),
     $  ern(nx),dlt4(nx),zr(nx),vr(nx),zs(nx),vs(nx),pmlts(nx),
     $  pmltg(nx),zg(nx),vg(nx),qsi(nx),ssi(nx),esi(nx),esw(nx),qsw(nx),
     $  ssw(nx),pihom(nx),pimlt(nx),dda(nb)
      common/bsat1/ pidw(nx),psaut(nx),qracs(nx),
     $  psaci(nx),psacw(nx),qsacw(nx),praci(nx),piacr(nx),praut(nx),
     $  pracw(nx),psfw(nx),psfi(nx),dgacs(nx),dgacw(nx),dgaci(nx),
     $  dgacr(nx),pgacs(nx),wgacs(nx),qgacw(nx),wgaci(nx),qgacr(nx),
     $  pgwet(nx),pgaut(nx),pracs(nx),psacr(nx),qsacr(nx),pgfr(nx),
     $  psmlt(nx),pgmlt(nx),psdep(nx),pgdep(nx),ddd(nb)
      common/badv1/ pt(nx),qv(nx),qc(nx),qr(nx),qi(nx),qs(nx),qg(nx),
     $        ddb(nb1)
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     $   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/brh1/ srro(nz),qrro(nz),sqc(nz),sqr(nz),sqi(nz),sqs(nz),
     $   sqg(nz),stqc(nz),stqr(nz),stqi(nz),stqs(nz),stqg(nz),ttt(nt5)
      common/ba/ y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),tair(nx),tairc(nx),
     $   pr(nx),ps(nx),pg(nx),prn(nx),psn(nx),dlt1(nx),dlt2(nx),
     $   dlt3(nx),rtair(nx)
      common/bsts/ thom(nz,4,7),tdw(nz,4,7),tmlt(nz,4,7),saut(nz,4,7),
     $ saci(nz,4,7),sacw(nz,4,7),raci(nz,4,7),tacr(nz,4,7),raut(nz,4,7),
     $ racw(nz,4,7),sfw(nz,4,7),sfi(nz,4,7),gacs(nz,4,7),gacw(nz,4,7),
     $ gaci(nz,4,7),gacr(nz,4,7),gwet(nz,4,7),gaut(nz,4,7),racs(nz,4,7),
     $ sacr(nz,4,7),gfr(nz,4,7),smlt(nz,4,7),gmlt(nz,4,7),sdep(nz,4,7),
     $ ssub(nz,4,7),gsub(nz,4,7),pern(nz,4,7),d3ri(nz,4,7),d3ir(nz,4,7),
     $ d2sr(nz,4,7),d2rs(nz,4,7),gdry(nz,4,7),coc(nz,4,7),coe(nz,4,7),
     $ smf0(nz,4,7),qc0(nz,4,7),qr0(nz,4,7),qi0(nz,4,7),qs0(nz,4,7),
     $ qg0(nz,4,7),sqc0(nz,4,7),sqr0(nz,4,7),sqi0(nz,4,7),sqs0(nz,4,7),
     $ sqg0(nz,4,7),erns(nz,4,7),wgrs(nz,4,7),qsws(nz,4,7),
     $ tbt(nz,4),qbt(nz,4)
      common/stls/ srsw(nz,4,7),srlw(nz,4,7)
      COMMON/BCS/ S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     $       S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     $       S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     $       SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     $       SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     $       SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCSS/ SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     $ SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     $ SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     $ SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     $ SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     $ SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)
      common/bch/ it(nx),iv(nt),ics(nx,4),ibz(nx,4)
      common/bch1/ rby(7)
      real y0(nx),zths(nx,itt),zqsTT(nx,itt),ztsTT(nx,itt),zpss(nx,itt)
      common/bls/ y0,zths,zqsTT,ztsTT,zpss
      real dbz(nx),fv(nz)
c-----------------------------------------------------------------------
cc    ***   three classes of ice-phase   (r&h)   **********************
       D22T=D2T
       IF(IJKADV .EQ. 0) THEN
         D2T=D2T
       ELSE
         D2T=DT
       ENDIF
        UCOR=3071.29/TNW**.75
      ft=dt/d2t
      rft=ril2*ft
      a0=ndt_stat
      do 10 i=1,imax
   10  it(i)=1
        r23t=rn23*d2t
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 1000 k=kt1,kt2
c        if (ijkadv .eq. 1) then
c          tb0=ta(k)
c          qb0=qa(k)
c        else
          tb0=ta1(k)
          qb0=qa1(k)
c        endif
c       p00=p0(k)
       rp0=3.799052e3/p0(k)
       pi0=pi(k)
       pir=1./(pi(k))
crh    pr0=1./p0(k)
       r00=rho(k)
       rr0=rrho(k)
       rrs=srro(k)
       rrq=qrro(k)
       fv0=fv(k)
       fvs=sqrt(fv(k))
          cp409=c409*pi0
          cv409=c409*avc
crh       alvr=r00*alv
          avcp=avc*pir
         zrr=1.e5*zrc*rrq
cs         vgcf=vgc*rrs
         r22f=rn22*fv0
         r231r=rn231*rr0
         r232rf=rn232*rrs*fvs
      do 125 i=it1,it2
       pt(i)=dpt(i,k)
       qv(i)=dqv(i,k)
       qc(i)=qcl(i,k)
       qr(i)=qrn(i,k)
c        if (qv(i)+qb0 .le. 0.) qv(i)=-qb0
        if (qc(i) .le. 1.e-20) qc(i)=0.0
        if (qr(i) .le. 1.e-20) qr(i)=0.0
       tair(i)=(pt(i)+tb0)*pi0
       tairc(i)=tair(i)-t0
        zr(i)=zrr
        vr(i)=0.0
        if (qr(i) .gt. 1.e-20) then
          dd(i)=r00*qr(i)
           y1(i)=sqrt(dd(i))
           y2(i)=sqrt(y1(i))
         zr(i)=zrc/y2(i)
         vr(i)=fv0*(vr0+vr1*y2(i)+vr2*y1(i)+vr3*y1(i)*y2(i))
          vr(i)=max(vr(i), 0.0)
        endif
         if (qr(i) .le. 1.e-12) vr(i)=0.0
c* 21 * praut   autoconversion of qc to qr                        **21**
c* 22 * pracw : accretion of qc by qr                             **22**
cso     praut(i)=max( rd1*(qc(i)-bound) , 0.0E0)
cso      pracw(i)=0.0
cso     if (qr(i) .gt. 1.e-16) then
cso      pracw(i)=rd2*qc(i)*qr(i)**.875
cso     endif
        praut(i)=max(rn21*(qc(i)-bnd21),0.0E0)
          y1(i)=1./zr(i)
         y2(i)=y1(i)*y1(i)
         y3(i)=y1(i)*y2(i)
         y4(i)=r22f*qc(i)*y3(i)*(rn50+rn51*y1(i)+rn52*y2(i)+rn53*y3(i))
        pracw(i)=max(y4(i),0.0E0)
c********   handling the negative cloud water (qc)    ******************
        y1(i)=d2t*(praut(i)+pracw(i))
       qc(i)=qc(i)-y1(i)
        y3(i)=qc(i)
        if (y3(i) .lt. 0.0) then
           y2(i)=1.
            if (y1(i) .ne. 0.0) y2(i)=qc(i)/y1(i)+1.
          praut(i)=praut(i)*y2(i)
          pracw(i)=pracw(i)*y2(i)
          qc(i)=0.0
         endif
        qr(i)=qr(i)+d2t*(praut(i)+pracw(i))
  125 continue
c*****   tao et al (1989) saturation technique  ***********************
        do 275 i=it1,it2
          tair(i)=(pt(i)+tb0)*pi0
          Y1(I)=1./(TAIR(I)-C358)
          QSW(I)=RP0*EXP(C172-C409*Y1(I))
          DD(I)=CP409*Y1(I)*Y1(I)
          DM(I)=QV(I)+QB0-QSW(I)
          CND(I)=DM(I)/(1.+AVCP*DD(I)*QSW(I))
CC    ******   CONDENSATION OR EVAPORATION OF QC  ******
          CND(I)=MAX(-QC(I), CND(I))
         PT(I)=PT(I)+AVCP*CND(I)
         QV(I)=QV(I)-CND(I)
         QC(I)=QC(I)+CND(I)
c* 23 * ern : evaporation of qr                                   **23**
        ern(i)=0.0
        if (qr(i) .gt. 1.e-20) then
         tair(i)=(pt(i)+tb0)*pi0
          rtair(i)=1./(tair(i)-c358)
           y2(i)=exp( c172-c409*rtair(i) )
          esw(i)=c610*y2(i)
          qsw(i)=rp0*y2(i)
          ssw(i)=(qv(i)+qb0)/qsw(i)-1.
          dm(i)=qv(i)+qb0-qsw(i)
           rsub1(i)=cv409*qsw(i)*rtair(i)*rtair(i)
          dd1(i)=max(-dm(i)/(1.+rsub1(i)),0.0E0)
cso         y1(i)=r00*qrn(i,k)
cso       ern(i)=(((1.6+124.9*y1(i)**.2046)*y1(i)**.525)
cso  1          /(2.55e6/(p00*qsw(i))+5.4e5))*(-dm(i)/(r00*qsw(i)))*d2t
            y3(i)=1./tair(i)
           dd(i)=y3(i)*(rn30a*y3(i)-rn10b)+rn10c*tair(i)/esw(i)
           y1(i)=-r23t*ssw(i)*(r231r/zr(i)**2+r232rf/zr(i)**3)/dd(i)
          ern(i)=min(dd1(i),qr(i),max(y1(i),0.0E0))
          pt(i)=pt(i)-avcp*ern(i)
          qv(i)=qv(i)+ern(i)
          qr(i)=qr(i)-ern(i)
         endif
c        if (qv(i)+qb0 .le. 0.) qv(i)=-qb0
        if(qc(i).lt.1.e-20) qc(i)=0.0
        if(qr(i).lt.1.e-20) qr(i)=0.0
       dpt(i,k)=pt(i)
       dqv(i,k)=qv(i)
       qcl(i,k)=qc(i)
       qrn(i,k)=qr(i)
  275  continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       scc=0.
       see=0.
       do 105 i=it1,it2
         dd(i)=max(-cnd(i), 0.e0)
         cnd(i)=max(cnd(i), 0.e0)
         scc=scc+cnd(i)
  105    see=see+dd(i)+ern(i)
        sc(k)=scc+sc(k)
        se(k)=see+se(k)
cc    ***   statistics for convective and anvil regimes   ***********
       if (id .eq. 1) then
        do 110 i=it1,it2
         dda(i)=rsw(i,k)
         ddb(i)=rlw(i,k)
         dm(i)=a0*(rho1(k)*ww1(i,k)+rho1(k+1)*ww1(i,k+1)+
     1              Y0(I)*(rho1(k)*wb(k)+rho1(k+1)*wb(k+1)))
  110    RQ(I)=.005*(RHO1(K)*(WW1(I,K)+wb(k))+
     1               RHO1(K+1)*(WW1(I,K+1)+wb(k+1)))/r00

        do 1500 kc=1,7
         do 320 mt=1,4
         do 320 i=it1,it2
          ibz(i,mt)=0
  320      if(ics(i,mt).eq.1) ibz(i,mt)=1
         do 325 i=it1,it2
          y1(i)=qc(i)+qr(i)
cc        y1(i)=qr(i)+qs(i)+qg(i)
cc         if(y1(i).ge.1.e-8 .and. ics(i,2).ne.1) then
c          if(y1(i).ge.1.e-20) then
            ibz(i,1)=1
c          endif
  325     continue
         do 330 mt=1,4
         do 330 i=it1,it2
           if(kc.eq.4) go to 330
           if(kc.le.3) go to 36
            if (rq(i).gt.rby(kc)) ibz(i,mt)=0
            go to 330
   36       if (rq(i).lt.rby(kc)) ibz(i,mt)=0
  330     continue
         do 350 mt=1,4
          sww=0.0
          scc=0.0
          see=0.0
          a1=0.0
          a9=0.0
          a10=0.0
          a11=0.0
          a25=0.0
          a26=0.0
          a27=0.0
          a33=0.0
          a34=0.0
          a38=0.0
          a41=0.0
          a42=0.0
          do 30 i=it1,it2
          if(ibz(i,mt).eq.1) then
          sww=sww+dm(i)
          scc=scc+cnd(i)
          see=see+dd(i)+ern(i)
          a9=a9+praut(i)
          a10=a10+pracw(i)
          a25=a25+dd(i)
          a26=a26+pgdep(i)
          a27=a27+dd1(i)
          a33=a33+qcl(i,k)
          a34=a34+qrn(i,k)
          a38=a38+ern(i)
          a41=a41+dda(i)
          a42=a42+ddb(i)
          endif
   30   continue
        smf0(k,mt,kc)=sww+smf0(k,mt,kc)
        coc(k,mt,kc)=scc+coc(k,mt,kc)
        coe(k,mt,kc)=see+coe(k,mt,kc)
        raut(k,mt,kc)=raut(k,mt,kc)+a9
        racw(k,mt,kc)=racw(k,mt,kc)+a10
        ssub(k,mt,kc)=ssub(k,mt,kc)+a25
        gsub(k,mt,kc)=gsub(k,mt,kc)+a26
        pern(k,mt,kc)=pern(k,mt,kc)+a27
        srsw(k,mt,kc)=srsw(k,mt,kc)+a41
        srlw(k,mt,kc)=srlw(k,mt,kc)+a42
        erns(k,mt,kc)=erns(k,mt,kc)+a38
        qc0(k,mt,kc)=a33*ril2
        qr0(k,mt,kc)=a34*ril2
        sqc0(k,mt,kc)=sqc0(k,mt,kc)+qc0(k,mt,kc)
        sqr0(k,mt,kc)=sqr0(k,mt,kc)+qr0(k,mt,kc)
  350   continue
 1500  continue
cc    ********************************************
c      do 40 i=2,iles
c       cnd(i)=cnd(i)+dep(i)
c       ern(i)=ern(i)+dd(i)+dd1(i)
c       cnd(i)=cnd(i)*rft
c       ern(i)=ern(i)*rft
c       ww=rq(i)
c        iww=max(ww,-5.E0)
c        iww=min(ww,15.E0)
c       if(ww.ge.0.) scu1(k)=scu1(k)+cnd(i)
c       if(ww.lt.0.) sed1(k)=sed1(k)+ern(i)
c       if(iww.lt.1) go to 45
c        s9(iww,k)=s9(iww,k)+cnd(i)
c        s10(iww,k)=s10(iww,k)+ern(i)
c       go to 40
c  45   jww=iabs(iww)
c       if(jww.lt.1) go to 40
c        sn9(jww,k)=sn9(jww,k)+cnd(i)
c        sn10(jww,k)=sn10(jww,k)+ern(i)
c  40   continue

C--------------------------------------------------------------------
C    CONDENSATION:  CND(IJ)
C    EVAPORATION:   DD(IJ)+ERN(IJ)
C    MASS FLUX:     DM(IJ)
C    CLOUD WATER:   QC(IJ)
C    RAIN:
C----------------------------------------------------------------------
       DO 42 IJ=it1,it2
        CND(IJ)=CND(IJ)*RFT
        ERN(IJ)=(ERN(IJ)+DD(IJ))*RFT
        Y5(IJ)=DM(IJ)
        A1=1.E6*R00*QR(IJ)

        A11=UCOR*(MAX(1.E-5,A1))**1.75
        ZDRY=MAX(1.,A11)
        DBZ(IJ)=10.*LOG10(ZDRY)
        DBZ(IJ)=DBZ(IJ)
        QC(IJ)=QC(IJ)
        QR(IJ)=QR(IJ)
   42  CONTINUE
        DO 44 IJ=it1,it2
        IF(RQ(IJ) .GE. 0.) SCU1(K)=SCU1(K)+CND(IJ)
        IF(RQ(IJ) .LT. 0.) SED1(K)=SED1(K)+ERN(IJ)
   44  CONTINUE
       DO 40 IJ=it1,it2
          IWW=0
          JWW=0
          IF (RQ(IJ).GT.-0.5) IWW=MIN(RQ(IJ)+0.5,15.)+1
          IF (RQ(IJ).LE.-0.5) JWW=MAX(RQ(IJ)-0.5,-5.)
          JWW=IABS(JWW)
        IF (ICS(IJ,2) .EQ. 1 .AND. IWW .GE. 1) THEN
          S9(IWW,K)=S9(IWW,K)+CND(IJ)
          S10(IWW,K)=S10(IWW,K)+ERN(IJ)
          S15(IWW,K)=S15(IWW,K)+Y5(IJ)
          S16(IWW,K)=S16(IWW,K)+QC(IJ)
          S17(IWW,K)=S17(IWW,K)+QR(IJ)
          S21(IWW,K)=S21(IWW,K)+DBZ(IJ)
          SCNT(IWW,K)=SCNT(IWW,K)+1.
        ENDIF
        IF (ICS(IJ,2) .EQ. 1 .AND. JWW .GE. 1) THEN
          SN9(jww,K)=SN9(jww,K)+CND(IJ)
          SN10(jww,K)=SN10(jww,K)+ERN(IJ)
          SN15(jww,K)=SN15(jww,K)+Y5(IJ)
          SN16(jww,K)=SN16(jww,K)+QC(IJ)
          SN17(jww,K)=SN17(jww,K)+QR(IJ)
          SN21(jww,K)=SN21(jww,K)+DBZ(IJ)
          SNCNT(jww,K)=SNCNT(jww,K)+1.
        ENDIF
        IF (ICS(IJ,3) .EQ. 1 .AND. IWW .GE. IWW) THEN
          SS9(IWW,K)=SS9(IWW,K)+CND(IJ)
          SS10(IWW,K)=SS10(IWW,K)+ERN(IJ)
          SS15(IWW,K)=SS15(IWW,K)+Y5(IJ)
          SS16(IWW,K)=SS16(IWW,K)+QC(IJ)
          SS17(IWW,K)=SS17(IWW,K)+QR(IJ)
          SS21(IWW,K)=SS21(IWW,K)+DBZ(IJ)
          SSCNT(IWW,K)=SSCNT(IWW,K)+1.
        ENDIF
        IF(ICS(IJ,3) .EQ. 1 .AND. jww.GE.1) THEN
          SSN9(jww,K)=SSN9(jww,K)+CND(IJ)
          SSN10(jww,K)=SSN10(jww,K)+ERN(IJ)
          SSN15(jww,K)=SSN15(jww,K)+Y5(IJ)
          SSN16(jww,K)=SSN16(jww,K)+QC(IJ)
          SSN17(jww,K)=SSN17(jww,K)+QR(IJ)
          SSN21(jww,K)=SSN21(jww,K)+DBZ(IJ)
          SSNCNT(jww,K)=SSNCNT(jww,K)+1.
        ENDIF
   40  CONTINUE

       ENDIF
 1000 continue
c     ****************************************************************
      if (id .eq. 1) then
       lc=0
       ls=0
       ln=0
       do 390 i=2,iles
        lc=lc+ics(i,2)
        ls=ls+ics(i,3)
        ln=ln+ics(i,4)
  390  continue
        if (lc .eq. 0) lc=1000000
        if (ls .eq. 0) ls=1000000
        if (ln .eq. 0) ln=1000000
       do 400 mt=1,4
        a1=ril2
         if (mt .eq. 2) a1=1./float(lc)
         if (mt .eq. 3) a1=1./float(ls)
         if (mt .eq. 4) a1=1./float(ln)
       do 410 k=2,kles
        y1(k)=0.0
        y2(k)=0.0
       do 410 i=2,iles
        if(ics(i,mt).eq.1) then
         y1(k)=y1(k)+dpt(i,k)
         y2(k)=y2(k)+dqv(i,k)
        endif
  410  continue
       do 430 k=2,kles
         tbt(k,mt)=y1(k)*a1
  430    qbt(k,mt)=y2(k)*a1
  400 continue
      endif
c     ****************************************************************
      do 460 k=2,kles
       y1(k)=0.
       y2(k)=0.
       sq(k)=0.
       sqc(k)=0.
       sqr(k)=0.
       sqi(k)=0.
       sqs(k)=0.
  460  sqg(k)=0.
      do 465 i=2,iles
      do 465 k=2,kles
       y1(k)=y1(k)+dpt(i,k)
       y2(k)=y2(k)+dqv(i,k)
       sq(k)=sq(k)+qcl(i,k)+qrn(i,k)
       sqc(k)=sqc(k)+qcl(i,k)
       sqr(k)=sqr(k)+qrn(i,k)
  465 continue
      do 470 k=2,kles
       y1(k)=y1(k)*ril2
       y2(k)=y2(k)*ril2
       sq(k)=sq(k)*ril2
       sqc(k)=sqc(k)*ril2
       sqr(k)=sqr(k)*ril2
       stqc(k)=stqc(k)+sqc(k)
       stqr(k)=stqr(k)+sqr(k)
  470 continue
      do 500 k=2,kles
       do 510 i=2,iles
        dpt(i,k)=dpt(i,k)-y1(k)
        dqv(i,k)=dqv(i,k)-y2(k)
  510  continue
       if(n.gt.1) go to 515
         IF (IJKADV .EQ. 1) THEN
           Y3(K)=TA(K)
           Y4(K)=QA(K)
          TA(K)=TA(K)+Y1(K)
          QA(K)=QA(K)+Y2(K)
          TA1(K)=y3(K)
          QA1(K)=y4(K)
         ELSE
          TA(K)=TA1(K)+Y1(K)
          QA(K)=QA1(K)+Y2(K)
         ENDIF
       go to 516
  515  do 517 i=2,iles
        IF (IJKADV .EQ. 0) THEN
         DPT1(I,K)=DPT1(I,K)+EPS*DPT(I,K)
         DQV1(I,K)=DQV1(I,K)+EPS*DQV(I,K)
         QCL1(I,K)=QCL1(I,K)+EPS*QCL(I,K)
         QRN1(I,K)=QRN1(I,K)+EPS*QRN(I,K)
c        ELSE
c          DPT1(I,K)=DPT(I,K)
c          DQV1(I,K)=DQV(I,K)
c          QCL1(I,K)=QCL(I,K)
c          QRN1(I,K)=QRN(I,K)
        ENDIF
  517  continue
         IF (IJKADV .EQ. 1) THEN
           Y3(K)=TA(K)
           Y4(K)=QA(K)
          TA(K)=TA(K)+Y1(K)
          QA(K)=QA(K)+Y2(K)
          TA1(K)=y3(K)
          QA1(K)=y4(K)
         ELSE
           A0=TA(K)+EPS*(-2.*TA(K)+TA1(K))
          TA(K)=TA1(K)+Y1(K)
          TA1(K)=A0+EPS*TA(K)
           A0=QA(K)+EPS*(-2.*QA(K)+QA1(K))
          QA(K)=QA1(K)+Y2(K)
          QA1(K)=A0+EPS*QA(K)
         ENDIF
  516  st(k)=ta(k)-tb(k)
       sv(k)=qa(k)-qb(k)
  500 continue
c     ****   fix the boundary conditions
      call boundy (dpt,dpt1)
      call boundy (dqv,dqv1)
      call boundy (qcl,qcl1)
      call boundy (qrn,qrn1)
c     ***********************************
      sq(1)=sq(2)
      ta(1)=ta(2)
      qa(1)=qa(2)
      ta1(1)=ta1(2)
      qa1(1)=qa1(2)
      sq(kmax)=sq(kles)
      ta(kmax)=ta(kles)*2.-ta(kl2)
      qa(kmax)=qa(kles)*2.-qa(kl2)
      ta1(kmax)=ta1(kles)*2.-ta1(kl2)
      qa1(kmax)=qa1(kles)*2.-qa1(kl2)

      do 600 k=3,kles
       IF (IJKADV .EQ. 1) THEN
        FD(K)=(TA1(K)-TA1(K-1))*RHO1(K)*.5
        FE(K)=(QA1(K)-QA1(K-1))*RHO1(K)*.5
       ELSE
        FD(K)=(TA(K)-TA(K-1))*RHO1(K)*.5
        FE(K)=(QA(K)-QA(K-1))*RHO1(K)*.5
       ENDIF
  600 continue
      fd(kmax)=fd(kles)
      fe(kmax)=fe(kles)
       d2t=d22t
      return
      end
