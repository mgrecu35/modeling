cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine satice (n,id,KCPS,fv)
c     ****   compute ice phase microphysics and saturation processes
      parameter (NX=514,NZ=43,NT=2880,ITT=244)
      parameter (nxi=nx-2,nzk=nz-2,nt5=5*nt)
      parameter (nb=nx*nz-31*nx,nb1=nx*nz-7*nx)
      common/timestat/ ndt_stat,itime_ave,mmave
      common/iice/ new_ice_sat
      common/itoga/ itoga,ISFC,ICE,ICE2
      common/iceopt/ ice913,ilif
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
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
      common/b1i/ qci(nx,nz)
      common/b1s/ qcs(nx,nz)
      common/b1g/ qcg(nx,nz)
      common/b2i/ qci1(nx,nz)
      common/b2s/ qcs1(nx,nz)
      common/b2g/ qcg1(nx,nz)
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
      common/b66b/ s_dep(nz),s_sub(nz),s_qrs(nz),s_qrl(nz),s_mel(nz),
     1   s_frz(nz)
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
c      common/bstsi/ ceds1i(nx,nz,4),ceds2i(nx,nz,4),
c     1   tstfi(nx,nz,4),tsqfi(nx,nz,4),rsli(nx,nz,4)
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
      real y0(nx),zths(nx,itt),zqsTT(nx,itt),ztsTT(nx,itt),zpss(nx,itt)
      common/bls/ y0,zths,zqsTT,ztsTT,zpss
      common/bch1/ rby(7)

c     for BM cumulus parameterization scheme
      common/cpsbm/ icps_bm,iexplicit
      COMMON/BPBL/ UHT(NZ),WHT(NZ),TGBAT0
      common/dinrad/ p00o(nz),dz0o(nz),tairsfco(nx),qairsfco(nx),
     1               thairsfco(nx),pairsfc(nx)
c
      common/sue/ ppress(nx,nz)
      common/sue1/ tmodo(nx,nz),qmodo(nx,nz),rainco(nx,1)
      common/sue3/ t_bm(nz,4,7),q_bm(nz,4,7)
C
c      COMMON/Q_BUGT/ Q1_G_H(NX,NZ),Q1_G_V(NX,NZ),
c     1               Q1A_G_H(NX,NZ),Q1A_G_V(NX,NZ),
c     2               Q1_D_H(NX,NZ),Q1_D_V(NX,NZ),
c     3               Q1A_D_H(NX,NZ),Q1A_D_V(NX,NZ),
c     4               Q2_G_H(NX,NZ),Q2_G_V(NX,NZ),
c     5               Q2A_G_H(NX,NZ),Q2A_G_V(NX,NZ),
c     6               Q2_D_H(NX,NZ),Q2_D_V(NX,NZ),
c     7               Q2A_D_H(NX,NZ),Q2A_D_V(NX,NZ),
c     8               Q1_HYD(NX,NZ),Q2_HYD(NX,NZ),
c     9               Q1A_HYD(NX,NZ),Q2A_HYD(NX,NZ),
c     9               Q1_RAD(NX,NZ),Q1A_RAD(NX,NZ),
c     9               IBUDSEC,RBUD
C
c
      common /lhblock/ rlh(nx,nz)
      dimension sm(nxi,1)
cc      dimension xland(nx,1),tref(nx,nzk)
      dimension plow(nxi,1),ptop(nxi,1),rainc(nxi,1)
      dimension zzh(nxi,1,nzk),tmod(nxi,nzk),qmod(nxi,nzk)
      dimension temp(nxi,1,nzk),xmoist(nxi,1,nzk)
      dimension preb(nxi,1,nzk),preh(nxi,1,nzk)
c
c      real tttbud(nx),qqqbud(nx)
      real    aa1(31),aa2(31),dbz(nx),fv(nz),cnd1(nx),dep1(nx)
         data aa1/.7939e-7,.7841e-6,.3369e-5,.4336e-5,.5285e-5,.3728e-5,
     1   .1852e-5,.2991e-6,.4248e-6,.7434e-6,.1812e-5,.4394e-5,.9145e-5,
     2   .1725e-4,.3348e-4,.1725e-4,.9175e-5,.4412e-5,.2252e-5,.9115e-6,
     3   .4876e-6,.3473e-6,.4758e-6,.6306e-6,.8573e-6,.7868e-6,.7192e-6,
     4   .6513e-6,.5956e-6,.5333e-6,.4834e-6/
        data aa2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
c-----------------------------------------------------------------------
      save
cc    ***   three classes of ice-phase   (r&h)   **********************
       D22T=D2T
       IF(IJKADV .EQ. 0) THEN
         D2T=D2T
       ELSE
         D2T=DT
       ENDIF
c
       call sat_zero
c
       cmin=1.e-20
       cmin1=1.e-16
c       cmin2=1.e-30
c       cmax2=1.e20
c
       R2ICE=1.
       IF (ICE2 .EQ. 1) R2ICE=0.
        UCOR=3071.29/TNW**.75                                            
        UCOS=687.97*ROQS**.25/TNS**.75
        UCOG=687.97*ROQG**.25/TNG**.75
        UWET=4.464**.95
      rt0=1./(t0-t00)
      ft=dt/d2t
      rft=ril2*ft
      A0=ndt_stat*.5*RIL2
      do 10 i=1,imax
       cnd1(i)=0.
       dep1(i)=0.
   10  it(i)=1
        bs3=bs+3.
        bg3=bg+3.
        bsh5=2.5+bsh
        bgh5=2.5+bgh
        bs6=6.+bs
        betah=.5*beta
        rdt=1./d2t
        r10t=rn10*d2t
        r11t=rn11*d2t
        r19t=rn19*d2t
        r19at=rn19a*d2t
        r20t=rn20*d2t
        r23t=rn23*d2t
c         r25a=rn25  	!tao new

         r25a=rn25

        r30t=rn30*d2t
        r33t=rn33*d2t
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NEW with BM CPS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (icps_bm .EQ. 1) THEN
        CPRA=1.004/.287
c         do ii=1,nx
c            xland(ii,1)=0.0
c         enddo

c         pre_ref=ppress(2,kles)

         do i=2,iles
            im=i-1
           do k=kles,2,-1
              kmm=kmax-k
             temp(im,1,kmm)=(ta1(k)+dpt(i,k))*pi(k)
             xmoist(im,1,kmm)=qa1(k)+dqv(i,k)
             preb(im,1,kmm)=100.e3*exp(CPRA*log(pi(k)))
c             preh(im,1,kmm)=100.e3*exp(cpra*log((ppress(i,k)-pre_ref)
c     1                                 +pi(k)))
             preh(im,1,kmm)=preb(im,1,kmm)
             zzh(im,1,kmm)=UHT(K)
           enddo
         enddo
c
c...... sm(nx,1) is surface mask:  = 1 water;  = 0 land
c
         do i=2,iles
             im=i-1
           SM(im,1)=1.
           rainc(im,1)=0.
           ptop(im,1)=.1*p0(kmax)
           plow(im,1)=.1*pairsfc(i)
         enddo

         delt2=d2t
c
c        time scale in BM is 3000 s, we could modify it by .5 or .25
c
         tscale=1.0
c
         jslab=1
         ibeg=1
         iend=il2

c         CALL BMPARA (DELT2,KCPS,TEMP,XMOIST,TMOD,QMOD,RAINC,JSLAB,IBEG,
c     1                IEND,SM,PTOP,PLOW,PREH,PREB,ZZH)

         kcps=kcps+1
c
         do i=2,iles
              im=i-1
            rainco(i,1)=1000.*rainc(im,1)/delt2
         enddo
c
           RR_AVC=1./AVC
c
         do i=2,iles
              im=i-1
            do k = kles,2,-1
               kmm=kmax-k
              dpt(i,k)=dpt(i,k)+tscale*tmod(im,kmm)*delt2/pi(k)
              dqv(i,k)=dqv(i,k)+tscale*qmod(im,kmm)*delt2
c
              TMODO(I,K)=tscale*TMOD(IM,KMM)*DELT2*RR_AVC
              QMODO(I,K)=-tscale*QMOD(IM,KMM)*DELT2
c
            enddo
         enddo
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 1000 k=kt1,kt2
c
       IF (IEXPLICIT .EQ. 1) THEN
c
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
c       r0s=sqrt(rho(k))
       rr0=rrho(k)
       rrs=srro(k)
       rrq=qrro(k)
       fv0=fv(k)
       fvs=sqrt(fv(k))
          cp409=c409*pi0
          cv409=c409*avc
          cp580=c580*pi0
          cs580=c580*asc
crh       alvr=r00*alv
          afcp=afc*pir
          avcp=avc*pir
          ascp=asc*pir
         zrr=1.e5*zrc*rrq
         zsr=1.e5*zsc*rrq
         zgr=1.e5*zgc*rrq
         vscf=vsc*fv0
         vgcf=vgc*fv0
cs         vgcf=vgc*rrs
         r1r=rn1*rr0
         r3f=rn3*fv0
         r4f=rn4*fv0
         r5f=rn5*fv0
         r6f=rn6*fv0
         r7rf=rn7*rr0*fv0
         r8rf=rn8*rr0*fv0
         r9rf=rn9*rr0*fv0
         r101r=rn101*rr0
         r102rf=rn102*rrs*fvs
         r12r=rn12*r00
         r14f=rn14*fv0
         r15f=rn15*fv0
c         r15af=rn15a*fv0
         r16rf=rn16*rr0*fv0
c         r172r=rn172*rr0
cc         r17arf=rn17a*r0s*fvs
         r18r=rn18*rr0
         r191r=rn191*rr0
         r192rf=rn192*rrs*fvs
         r22f=rn22*fv0
         r231r=rn231*rr0
         r232rf=rn232*rrs*fvs
         r25rt=rn25*rr0*d2t
         r31r=rn31*rr0
         r32rt=rn32*d2t*rrs
         r331r=rn331*rr0
         r332rf=rn332*rrs*fvs
         r34f=rn34*fv0
      do 125 i=it1,it2

c       tttbud(i)=dpt(i,k)
c       qqqbud(i)=dqv(i,k)

       pt(i)=dpt(i,k)
       qv(i)=dqv(i,k)
       qc(i)=qcl(i,k)
       qr(i)=qrn(i,k)
       qi(i)=qci(i,k)
       qs(i)=qcs(i,k)
       qg(i)=qcg(i,k)
c        if (qv(i)+qb0 .le. 0.) qv(i)=-qb0
        if (qc(i) .le. CMIN) qc(i)=0.0
        if (qr(i) .le. CMIN) qr(i)=0.0
        if (qi(i) .le. CMIN) qi(i)=0.0
        if (qs(i) .le. CMIN) qs(i)=0.0
        if (qg(i) .le. CMIN) qg(i)=0.0
       tair(i)=(pt(i)+tb0)*pi0
       tairc(i)=tair(i)-t0
        zr(i)=zrr
        vr(i)=0.0
        zs(i)=zsr
        vs(i)=0.0
        zg(i)=zgr
        vg(i)=0.0
        if (qr(i) .gt. CMIN) then
          dd(i)=r00*qr(i)
           y1(i)=sqrt(dd(i))
           y2(i)=sqrt(y1(i))
         zr(i)=zrc/y2(i)
         vr(i)=fv0*(vr0+vr1*y2(i)+vr2*y1(i)+vr3*y1(i)*y2(i))
          vr(i)=max(vr(i), 0.0)
        endif
        if (qs(i) .gt. CMIN) then
          dd(i)=r00*qs(i)
           y1(i)=dd(i)**.25
         zs(i)=zsc/y1(i)
         vs(i)=max(vscf*dd(i)**bsq, 0.0)
        endif
        if (qg(i) .gt. CMIN) then
          dd(i)=r00*qg(i)
           y1(i)=dd(i)**.25
         zg(i)=zgc/y1(i)
         vg(i)=max(vgcf*dd(i)**bgq, 0.0)
        endif
          if (qr(i) .le. cmin1) vr(i)=0.0
          if (qs(i) .le. cmin1) vs(i)=0.0
          if (qg(i) .le. cmin1) vg(i)=0.0
  125   continue
c*  1 * psaut : autoconversion of qi to qs                        ***1**
c*  3 * psaci : accretion of qi to qs                             ***3**
c*  4 * psacw : accretion of qc by qs (riming) (qsacw for psmlt)  ***4**
c* 34 * pwacs : collection of qs by qc                            **34**
c*  5 * praci : accretion of qi by qr                             ***5**
c*  6 * piacr : accretion of qr or qg by qi                       ***6**
      do 150 i=it1,it2
       psaut(i)=0.0
       psaci(i)=0.0
       praci(i)=0.0
       piacr(i)=0.0
       psacw(i)=0.0
       pwacs(i)=0.0
       qsacw(i)=0.0
       if(tair(i).lt.t0) then
         y1(i)=rdt*(qi(i)-r1r*exp(beta*tairc(i)))
        psaut(i)=max(y1(i), 0.0E0)

        psaci(i)=r3f*qi(i)/zs(i)**bs3
        psacw(i)=r4f*qc(i)/zs(i)**bs3
        pwacs(i)=R2ICE*r34f*qc(i)/zs(i)**bs6
          y1(i)=1./zr(i)
         y2(i)=y1(i)*y1(i)
         y3(i)=y1(i)*y2(i)
         dd(i)=r5f*qi(i)*y3(i)*(rn50+rn51*y1(i)+rn52*y2(i)+rn53*y3(i))
        praci(i)=max(dd(i),0.0E0)
           y4(i)=y3(i)*y3(i)
         dd1(i)=r6f*qi(i)*y4(i)*(rn60+rn61*y1(i)+rn62*y2(i)+rn63*y3(i))
        piacr(i)=max(dd1(i),0.0E0)
       else
        qsacw(i)=r4f*qc(i)/zs(i)**bs3
       endif
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
c* 12 * psfw : bergeron processes for qs (koening, 1971)          **12**
c* 13 * psfi : bergeron processes for qs                          **13**
          psfw(i)=0.0
          psfi(i)=0.0
           if(tair(i).lt.t0) then
             y1(i)=max( min(tairc(i), -1.e0), -31.e0)
             it(i)=int(abs(y1(i)))
             y1(i)=rn12a(it(i))
             y2(i)=rn12b(it(i))
             y3(i)=rn13(it(i))
            psfw(i)=max(d2t*y1(i)*(y2(i)+r12r*qc(i))*qi(i),0.0e0)
            psfi(i)=y3(i)*qi(i)
           endif
cttt***** qg=qg+min(pgdry,pgwet)
c*  9 * pgacs : accretion of qs by qg (dgacs,wgacs: dry and wet)  ***9**
c* 14 * dgacw : accretion of qc by qg (qgacw for pgmlt)           **14**
c* 15 * dgaci : accretion of qi by qg (wgaci for wet growth)      **15**
c* 16 * dgacr : accretion of qr to qg (qgacr for pgmlt)           **16**
         y1(i)=abs( vg(i)-vs(i) )
         y2(i)=zs(i)*zg(i)
         y3(i)=5./y2(i)
         y4(i)=.08*y3(i)*y3(i)
         y5(i)=.05*y3(i)*y4(i)
         y2(i)=y1(i)*(y3(i)/zs(i)**5+y4(i)/zs(i)**3+y5(i)/zs(i))
        pgacs(i)=R2ICE*r9rf*y2(i)
        dgacs(i)=pgacs(i)
crh     wgacs(i)=10.*r9rf*y2(i)
        wgacs(i)=0.0
         y1(i)=1./zg(i)**bg3
        dgacw(i)=R2ICE*r14f*qc(i)*y1(i)
        qgacw(i)=r14f*qc(i)*y1(i)
        dgaci(i)=R2ICE*r15f*qi(i)*y1(i)
crh     wgaci(i)=r15af*qi(i)*y1(i)
        wgaci(i)=0.0
         y1(i)=abs( vg(i)-vr(i) )
         y2(i)=zr(i)*zg(i)
         y3(i)=5./y2(i)
         y4(i)=.08*y3(i)*y3(i)
         y5(i)=.05*y3(i)*y4(i)
         dd(i)=r16rf*y1(i)*(y3(i)/zr(i)**5+y4(i)/zr(i)**3+y5(i)/zr(i))
        dgacr(i)=R2ICE*max(dd(i),0.0E0)
        qgacr(i)=dgacr(i)
         if (tair(i) .ge. t0) then
          dgacs(i)=0.0
crh       wgacs(i)=0.0
          dgacw(i)=0.0
          dgaci(i)=0.0
crh       wgaci(i)=0.0
          dgacr(i)=0.0
         else
          pgacs(i)=0.0
          qgacw(i)=0.0
          qgacr(i)=0.0
         endif
  150  continue
c*******pgdry : dgacw+dgaci+dgacr+dgacs                           ******
c* 17 * pgwet : wet growth of qg                                  **17**
crh       pgwet(i)=0.0
crh       if (tair(i) .lt. t0) then
crh         y1(i)=1./(alf+rn17c*tairc(i))
crh         y2(i)=rp0-(qv(i)+qb0)
crh         y3(i)=.78/zg(i)**2+r17arf/zg(i)**bgh5
crh         y4(i)=rn171*y2(i)-r172r*tairc(i)
crh         dd(i)=y1(i)*(y4(i)*y3(i)+(wgaci(i)
crh  1                               +wgacs(i))*(alf+rn17b*tairc(i)))
crh        pgwet(i)=max(dd(i), 0.0E0)
crh       endif
c******** shed process (wgacr=pgwet-dgacw-wgaci-wgacs)
crh     wgacr(i)=pgwet(i)-dgacw(i)-wgaci(i)-wgacs(i)
crh      y2(i)=dgacw(i)+dgaci(i)+dgacr(i)+dgacs(i)
crh      if (pgwet(i) .ge. y2(i)) then
crh       wgacr(i)=0.0
crh       wgaci(i)=0.0
crh       wgacs(i)=0.0
crh      else
crh       dgacr(i)=0.0
crh       dgaci(i)=0.0
crh       dgacs(i)=0.0
crh      endif
c*******pgdry : dgacw+dgaci+dgacr+dgacs                           ******
c* 15 * dgaci : accretion of qi by qg (wgaci for wet growth)      **15**
c* 17 * pgwet : wet growth of qg                                  **17**
c       do 160 i=it1,it2
c        pgwet(i)=0.0
c  160  continue
c********   handling the negative cloud water (qc)    ******************
c********   handling the negative cloud ice (qi)      ******************
       do 175 i=it1,it2
c
        pgwet(i)=0.0
c
         if (ice913 .eq. 0) then
            y1(i)=qc(i)/d2t
           psacw(i)=min(y1(i), psacw(i))
           praut(i)=min(y1(i), praut(i))
           pracw(i)=min(y1(i), pracw(i))
           psfw(i)= min(y1(i), psfw(i))
           dgacw(i)=min(y1(i), dgacw(i))
           qsacw(i)=min(y1(i), qsacw(i))
           qgacw(i)=min(y1(i), qgacw(i))
         endif
c
        y1(i)=d2t*(psacw(i)+praut(i)+pracw(i)+psfw(i)+dgacw(i)
     1            +qsacw(i)+qgacw(i))
        qc(i)=qc(i)-y1(i)
c
        if (qc(i) .lt. 0.0) then
          y2(i)=1.
           if(y1(i) .ne. 0.) y2(i)=qc(i)/y1(i)+1.
          psacw(i)=psacw(i)*y2(i)
          praut(i)=praut(i)*y2(i)
          pracw(i)=pracw(i)*y2(i)
          psfw(i)=psfw(i)*y2(i)
          dgacw(i)=dgacw(i)*y2(i)
          qsacw(i)=qsacw(i)*y2(i)
          qgacw(i)=qgacw(i)*y2(i)
          qc(i)=0.0
        endif
c
         if (ice913 .eq. 0) then
            y1(i)=qi(i)/d2t
           psaut(i)=min(y1(i), psaut(i))
           psaci(i)=min(y1(i), psaci(i))
           praci(i)=min(y1(i), praci(i))
           psfi(i)= min(y1(i), psfi(i))
           dgaci(i)=min(y1(i), dgaci(i))
           wgaci(i)=min(y1(i), wgaci(i))
         endif
c
        y1(i)=d2t*(psaut(i)+psaci(i)+praci(i)+psfi(i)+dgaci(i)+wgaci(i))
c
       qi(i)=qi(i)-y1(i)
       if (qi(i).lt.0.0) then
         y2(i)=1.
           if(y1(i) .ne. 0.) y2(i)=qi(i)/y1(i)+1.
           psaut(i)=psaut(i)*y2(i)
           psaci(i)=psaci(i)*y2(i)
           praci(i)=praci(i)*y2(i)
           psfi(i)=psfi(i)*y2(i)
           dgaci(i)=dgaci(i)*y2(i)
           wgaci(i)=wgaci(i)*y2(i)
           qi(i)=0.0
        endif
c
         wgacr(i)=qgacr(i)+qgacw(i)
        dlt3(i)=0.0
         if (qr(i) .lt. 1.e-4) dlt3(i)=1.
        dlt4(i)=1.
         if (qc(i) .gt. 5.e-4) dlt4(i)=0.0
         if (qs(i) .le. 1.e-4) dlt4(i)=1.
          if (tair(i) .ge. t0) then
           dlt3(i)=0.0
           dlt4(i)=0.0
          endif
          IF (ICE2 .EQ. 1) THEN
            DLT3(I)=1.0
            DLT4(I)=1.0
          ENDIF
        pr(i)=d2t*(qsacw(i)+praut(i)+pracw(i)+wgacr(i)-qgacr(i))
        ps(i)=d2t*(psaut(i)+psaci(i)+dlt4(i)*psacw(i)+psfw(i)+psfi(i)
     1             +dlt3(i)*praci(i))
        pg(i)=d2t*((1.-dlt3(i))*praci(i)+dgaci(i)+wgaci(i)+dgacw(i)
     1             +(1.-dlt4(i))*psacw(i))
  175 continue
c*  7 * pracs : accretion of qs by qr                             ***7**
c*  8 * psacr : accretion of qr by qs (qsacr for psmlt)           ***8**
      do 200 i=it1,it2
         y1(i)=abs( vr(i)-vs(i) )
         y2(i)=zr(i)*zs(i)
         y3(i)=5./y2(i)
         y4(i)=.08*y3(i)*y3(i)
         y5(i)=.05*y3(i)*y4(i)
           dd(i)=y1(i)*(y3(i)/zs(i)**5+y4(i)/zs(i)**3+y5(i)/zs(i))
        pracs(i)=r7rf*max(dd(i), 0.0E0)
        qracs(i)=min(d2t*pracs(i),qs(i))
        pracs(i)=R2ICE*pracs(i)
          dd(i)=y1(i)*(y3(i)/zr(i)**5+y4(i)/zr(i)**3+y5(i)/zr(i))
        psacr(i)=r8rf*max(dd(i), 0.0E0)
        qsacr(i)=psacr(i)
         if (tair(i) .ge. t0)then
          pgaut(i)=0.0
          pracs(i)=0.0
          psacr(i)=0.0
         else
          qsacr(i)=0.0
          qracs(i)=0.0
         endif
c*  2 * pgaut : autoconversion of qs to qg                        ***2**
c* 18 * pgfr : freezing of qr to qg                               **18**
        pgaut(i)=0.0
        pgfr(i)=0.0
         if (tair(i) .lt. t0 .and. tair(i) .ge. 248.16) then
crh         y1(i)=exp(.09*tairc(i))
crh         y2(i)=rn2*y1(i)*(qs(i)-bnd2)
crh       pgaut(i)=max(y2(i),0.0E0)
           y1(i)=exp(rn18a*(t0-tair(i)))
          pgfr(i)=R2ICE*max(r18r*(y1(i)-1.)/zr(i)**7., 0.0E0)
         endif
  200 continue
c********   handling the negative rain water (qr)    *******************
c********   handling the negative snow (qs)          *******************
      do 225 i=it1,it2
c
        if (ice913 .eq. 0) then
           y1(i)=qr(i)/d2t
          piacr(i)=min(y1(i), piacr(i))
          dgacr(i)=min(y1(i), dgacr(i))
          psacr(i)=min(y1(i), psacr(i))
          pgfr(i)= min(y1(i), pgfr(i))
        endif
c
         y1(i)=(piacr(i)+dgacr(i)+psacr(i)+pgfr(i))*d2t
        qr(i)=qr(i)+pr(i)+qracs(i)-y1(i)

         if (qr(i).lt.0.0) then
           y2(i)=1.
            if(y1(i) .ne. 0.) y2(i)=qr(i)/y1(i)+1.
          piacr(i)=piacr(i)*y2(i)
          dgacr(i)=dgacr(i)*y2(i)
crh       wgacr(i)=wgacr(i)*y2(i)
          pgfr(i)=pgfr(i)*y2(i)
          psacr(i)=psacr(i)*y2(i)
crh       qgacr(i)=qgacr(i)*y2(i)
crh       qsacr(i)=qsacr(i)*y2(i)
          qr(i)=0.0
        endif
c
        dlt2(i)=1.
         if (qr(i) .gt. 1.e-4) dlt2(i)=0.
         if (qs(i) .le. 1.e-4) dlt2(i)=1.
         if (tair(i) .ge. t0) dlt2(i)=0.
          IF (ICE2 .EQ. 1) THEN
            DLT2(I)=1.0
          ENDIF


        if (ice913 .eq. 0) then
           y1(i)=qs(i)/d2t
          pgacs(i)=min(y1(i), pgacs(i))
          dgacs(i)=min(y1(i), dgacs(i))
          wgacs(i)=min(y1(i), wgacs(i))
          pgaut(i)=min(y1(i), pgaut(i))
          pracs(i)=min(y1(i), pracs(i))
          pwacs(i)=min(y1(i), pwacs(i))
        endif
c
        prn(i)=d2t*((1.-dlt3(i))*piacr(i)+dgacr(i)+pgfr(i)+(1.-dlt2(i))
     1              *psacr(i))
        ps(i)=ps(i)+d2t*(dlt3(i)*piacr(i)+dlt2(i)*psacr(i))
         pracs(i)=(1.-dlt2(i))*pracs(i)
         pwacs(i)=(1.-dlt4(i))*pwacs(i)
c
      psn(i)=d2t*(pgacs(i)+dgacs(i)+wgacs(i)+pgaut(i)+pracs(i)+pwacs(i))

       qs(i)=qs(i)+ps(i)-qracs(i)-psn(i)

         if (qs(i) .lt. 0.0) then
           y2(i)=1.
            if(psn(i) .ne. 0.) y2(i)=qs(i)/psn(i)+1.
           pgacs(i)=pgacs(i)*y2(i)
           dgacs(i)=dgacs(i)*y2(i)
           wgacs(i)=wgacs(i)*y2(i)
           pgaut(i)=pgaut(i)*y2(i)
           pracs(i)=pracs(i)*y2(i)
           pwacs(i)=pwacs(i)*y2(i)
           qs(i)=0.0
         endif
c
c           if (qracs(i) .ge. qs(i)) qracs(i)=qs(i)
c            qr(i)=qr(i)+qracs(i)
c            qs(i)=qs(i)-qracs(i)
c
      psn(i)=d2t*(pgacs(i)+dgacs(i)+wgacs(i)+pgaut(i)+pracs(i)+pwacs(i))
       qg(i)=qg(i)+R2ICE*(pg(i)+prn(i)+psn(i))
        y1(i)=d2t*(psacw(i)+psfw(i)+dgacw(i)+piacr(i)+dgacr(i)
     1             +psacr(i)+pgfr(i))-qracs(i)
       pt(i)=pt(i)+afcp*y1(i)       !LH
       rlh(i,k)=rlh(i,k)+afcp*y1(i)   !LH
  225 continue
c* 11 * psmlt : melting of qs                                     **11**
c* 19 * pgmlt : melting of qg to qr                               **19**
        do 250 i=it1,it2
          psmlt(i)=0.0
          pgmlt(i)=0.0
        tair(i)=(pt(i)+tb0)*pi0
        if (tair(i) .ge. t0) then
          tairc(i)=tair(i)-t0
          dd(i)=r11t*tairc(i)*(r101r/zs(i)**2+r102rf/zs(i)**bsh5)
         psmlt(i)=min(qs(i),max(dd(i),0.0E0))
           y2(i)=r191r/zg(i)**2+r192rf/zg(i)**bgh5
          dd1(i)=tairc(i)*(r19t*y2(i)+r19at*(qgacw(i)+qgacr(i)))
         pgmlt(i)=R2ICE*min(qg(i),max(dd1(i), 0.0E0))
         pt(i)=pt(i)-afcp*(psmlt(i)+pgmlt(i))      !LH
         rlh(i,k)=rlh(i,k)-afcp*(psmlt(i)+pgmlt(i))  !LH
         qr(i)=qr(i)+psmlt(i)+pgmlt(i)
         qs(i)=qs(i)-psmlt(i)
         qg(i)=qg(i)-pgmlt(i)
        endif
c* 24 * pihom : homogeneous freezing of qc to qi (t < t00)        **24**
c* 25 * pidw : deposition growth of qc to qi ( t0 < t <= t00)     **25**
c* 26 * pimlt : melting of qi to qc (t >= t0)                     **26**
        if (qc(i).le.cmin) qc(i)=0.0
        if (qi(i).le.cmin) qi(i)=0.0
          tair(i)=(pt(i)+tb0)*pi0
         if(tair(i).le.t00) then
          pihom(i)=qc(i)
         else
          pihom(i)=0.0
         endif
         if(tair(i).ge.t0) then
          pimlt(i)=qi(i)
         else
          pimlt(i)=0.0
         endif
         pidw(i)=0.0
         if (tair(i).lt.t0 .and. tair(i).gt.t00) then
           tairc(i)=tair(i)-t0
            y1(i)=max( min(tairc(i), -1.e0), -31.e0)
            it(i)=int(abs(y1(i)))
c LFO
c	dd(i)=r25rt*rn25a(it(i))*exp(0.5*abs(tairc(i)))
c KFL formula (1993)
           y3(i)=aa2(it(i))
           y4(i)=exp(abs(0.5*tairc(i)))
           y5(i)=aa1(it(i))*(r00*qi(i)/(r25a*y4(i)))**y3(i)
           dd(i)=r25rt*y5(i)*y4(i)
          pidw(i)=min(qc(i),dd(i))
         endif
          y1(i)=pihom(i)-pimlt(i)+pidw(i)
        pt(i)=pt(i)+afcp*y1(i)        !LH
        rlh(i,k)=rlh(i,k)+afcp*y1(i)    !LH
        qc(i)=qc(i)-y1(i)
        qi(i)=qi(i)+y1(i)
c* 31 * pint  : initiation of qi                                  **31**
c* 32 * pidep : deposition of qi                                  **32**
        pint(i)=0.0
        if (ice913 .eq. 0) then
          tair(i)=(pt(i)+tb0)*pi0
          if (tair(i) .lt. t0) then
            if (qi(i) .le. cmin) qi(I)=0.
             tairc(i)=tair(i)-t0
            dd(i)=r31r*exp(beta*tairc(i))
             rtair(i)=1./(tair(i)-c76)
             y2(i)=exp(c218-c580*rtair(i))
            qsi(i)=rp0*y2(i)
             esi(i)=c610*y2(i)
            ssi(i)=(qv(i)+qb0)/qsi(i)-1.
              dm(i)=max( (qv(i)+qb0-qsi(i)), 0.e0)
              rsub1(i)=cs580*qsi(i)*rtair(i)*rtair(i)
            dep(i)=dm(i)/(1.+rsub1(i))
            pint(i)=max(min(dd(i), dm(i)), 0.)
             y1(i)=1./tair(i)
             y2(i)=exp(betah*tairc(i))
             y3(i)=sqrt(qi(i))
             dd(i)=y1(i)*(rn10a*y1(i)-rn10b)+rn10c*tair(i)/esi(i)
            pidep(i)=max(r32rt*ssi(i)*y2(i)*y3(i)/dd(i), 0.e0)
            pint(i)=pint(i)+pidep(i)
            pint(i)=min(pint(i), dep(i))
cc             if (pint(i) .le. cmin) pint(i)=0.
cc             if (qi(i) .le. cmin) qi(I)=0.
            pt(i)=pt(i)+ascp*pint(i)       !LH
            rlh(i,k)=rlh(i,k)+ascp*pint(i)   !LH
            qv(i)=qv(i)-pint(i)
            qi(i)=qi(i)+pint(i)
          endif
        endif
  250   continue
*****   tao et al (1989) saturation technique  ***********************
        IF (NEW_ICE_SAT .EQ. 0) THEN ! cccshie by tao 5/3/01
        do 275 i=it1,it2
         tair(i)=(pt(i)+tb0)*pi0
        cnd(i)=rt0*(tair(i)-t00)
        dep(i)=rt0*(t0-tair(i))
          y1(i)=1./(tair(i)-c358)
          y2(i)=1./(tair(i)-c76)
         qsw(i)=rp0*exp(c172-c409*y1(i))
         qsi(i)=rp0*exp(c218-c580*y2(i))
          dd(i)=cp409*y1(i)*y1(i)
          dd1(i)=cp580*y2(i)*y2(i)
         if (qc(i).le.cmin) qc(i)=cmin
         if (qi(i).le.cmin) qi(i)=cmin
         if (tair(i).ge.t0) then
          dep(i)=0.0
          cnd(i)=1.
          qi(i)=0.0
         endif
         if (tair(i).lt.t00) then
          cnd(i)=0.0
          dep(i)=1.
          qc(i)=0.0
         endif
          y5(i)=avcp*cnd(i)+ascp*dep(i)
c          if (qc(i) .ge. cmin .or. qi(i) .ge. cmin) then
           y1(i)=qc(i)*qsw(i)/(qc(i)+qi(i))
           y2(i)=qi(i)*qsi(i)/(qc(i)+qi(i))
          y4(i)=dd(i)*y1(i)+dd1(i)*y2(i)
         qvs(i)=y1(i)+y2(i)
         rsub1(i)=(qv(i)+qb0-qvs(i))/(1.+y4(i)*y5(i))
        cnd(i)=cnd(i)*rsub1(i)
        dep(i)=dep(i)*rsub1(i)
         if (qc(i).le.cmin) qc(i)=0.
         if (qi(i).le.cmin) qi(i)=0.
cc    ******   condensation or evaporation of qc  ******
         cnd(i)=max(-qc(i),cnd(i))
cc    ******   deposition or sublimation of qi    ******
         dep(i)=max(-qi(i),dep(i))
        pt(i)=pt(i)+avcp*cnd(i)+ascp*dep(i)   !LH
        rlh(i,k)=rlh(i,k)+avcp*cnd(i)+ascp*dep(i) !LH
        qv(i)=qv(i)-cnd(i)-dep(i)
        qc(i)=qc(i)+cnd(i)
        qi(i)=qi(i)+dep(i)
  275   continue
        ENDIF ! cccshie by tao 5/3/01
c
cc
c
        IF (NEW_ICE_SAT .EQ. 1) THEN
           DO I=IT1,IT2
              TAIR(I)=(PT(I)+TB0)*PI0
c             CND1(I)=RT0*(TAIR(I)-T00)
c             DEP1(I)=RT0*(T0-TAIR(I))
              CND(I)=RT0*(TAIR(I)-T00)  ! cccshie by tao 5/3/01
              DEP(I)=RT0*(T0-TAIR(I))   ! cccshie by tao 5/3/01
                Y1(I)=1./(TAIR(I)-C358)
                Y2(I)=1./(TAIR(I)-C76)
              QSW(I)=RP0*EXP(C172-C409*Y1(I))
              QSI(I)=RP0*EXP(C218-C580*Y2(I))
              DD(I)=CP409*Y1(I)*Y1(I)
              DD1(I)=CP580*Y2(I)*Y2(I)
              Y5(I)=AVCP*CND(I)+ASCP*DEP(I) ! cccshie by tao 5/3/01
              Y1(I)=RT0*(TAIR(I)-T00)*QSW(I) ! cccshie by tao 5/3/01
              Y2(I)=RT0*(T0-TAIR(I))*QSI(I)  ! cccshie by tao 5/3/01
c             IF (QC(I).LE.CMIN) QC(I)=CMIN  ! cccshie by tao 5/3/01
c             IF (QI(I).LE.CMIN) QI(I)=CMIN  ! cccshie by tao 5/3/01
              IF (TAIR(I).GE.T0) THEN
c                 DEP1(I)=0.0
c                 CND1(I)=1.
c                 QI(I)=0.0
                  DEP(I)=0.0 ! cccshie by tao 5/3/01
                  CND(I)=1.  ! cccshie by tao 5/3/01
                  Y2(I)=0.   ! cccshie by tao 5/3/01
                  Y1(I)=QSW(I) ! cccshie by tao 5/3/01
              ENDIF
              IF (TAIR(I).LT.T00) THEN
                 CND(I)=0.0  ! cccshie by tao 5/3/01
                 DEP(I)=1.   ! cccshie by tao 5/3/01
                 Y2(I)=QSI(I) ! cccshie by tao 5/3/01
                 Y1(I)=0.     ! cccshie by tao 5/3/01
c                CND1(I)=0.0
c                DEP1(I)=1.
c                QC(I)=0.0
              ENDIF
c             Y5(I)=AVCP*CND1(I)+ASCP*DEP1(I) ! cccshie by tao 5/3/01
c             Y1(I)=QC(I)*QSW(I)/(QC(I)+QI(I)) ! cccshie by tao 5/3/01
c             Y2(I)=QI(I)*QSI(I)/(QC(I)+QI(I)) ! cccshie by tao 5/3/01
              Y4(I)=DD(I)*Y1(I)+DD1(I)*Y2(I)
              QVS(I)=Y1(I)+Y2(I)
              RSUB1(I)=(QV(I)+QB0-QVS(I))/(1.+Y4(I)*Y5(I))
             CND(I)=CND(I)*RSUB1(I) ! cccshie by tao 5/3/01
             DEP(I)=DEP(I)*RSUB1(I) ! cccshie by tao 5/3/01
c            CND1(I)=CND1(I)*RSUB1(I)
c            DEP1(I)=DEP1(I)*RSUB1(I)
c             IF (QC(I).LE.CMIN) QC(I)=0. ! cccshie by tao 5/3/01
c             IF (QI(I).LE.CMIN) QI(I)=0. ! cccshie by tao 5/3/01
CC    ******   CONDENSATION OR EVAPORATION OF QC  ******
c            CND1(I)=MAX(-QC(I),CND1(I))
             CND(I)=MAX(-QC(I),CND(I)) ! cccshie by tao 5/3/01
CC    ******   DEPOSITION OR SUBLIMATION OF QI    ******
             DEP(I)=MAX(-QI(I),DEP(I)) ! cccshie by tao 5/3/01
             PT(I)=PT(I)+AVCP*CND(I)+ASCP*DEP(I) ! cccshie by tao 5/3/01     LH
             rlh(i,k)=rlh(i,k)+AVCP*CND(I)+ASCP*DEP(I) ! cccshie by tao 5/3/01 LH
             QV(I)=QV(I)-CND(I)-DEP(I) ! cccshie by tao 5/3/01
             QC(I)=QC(I)+CND(I) ! cccshie by tao 5/3/01
             QI(I)=QI(I)+DEP(I) ! cccshie by tao 5/3/01
c            DEP1(I)=MAX(-QI(I),DEP1(I))
c            PT(I)=PT(I)+AVCP*CND1(I)+ASCP*DEP1(I)
c            QV(I)=QV(I)-CND1(I)-DEP1(I)
c            QC(I)=QC(I)+CND1(I)
c            QI(I)=QI(I)+DEP1(I)
           ENDDO
        ENDIF
c
cc
c
c* 10 * psdep : deposition of qs                                  **10**
c* 20 * pgdep : deposition of qg                                  **20**
        tpsg=1.e-8
          if (ice913 .eq. 0) tpsg=1.e-5
        do 300 i=it1,it2
       psdep(i)=0.0
       pgdep(i)=0.0
       tair(i)=(pt(i)+tb0)*pi0
        if (tair(i) .lt. t0) then
          if(qc(i)+qi(i) .gt. tpsg) then
           dlt1(i)=1.
          else
           dlt1(i)=0.
          endif
         rtair(i)=1./(tair(i)-c76)
          y2(i)=exp(c218-c580*rtair(i))
         qsi(i)=rp0*y2(i)
         esi(i)=c610*y2(i)
         ssi(i)=dlt1(i)*((qv(i)+qb0)/qsi(i)-1.)
          dm(i)=qv(i)+qb0-qsi(i)
          rsub1(i)=cs580*qsi(i)*rtair(i)*rtair(i)
          dd1(i)=max(dm(i)/(1.+rsub1(i)),0.0e0)
           y3(i)=1./tair(i)
          dd(i)=y3(i)*(rn10a*y3(i)-rn10b)+rn10c*tair(i)/esi(i)
           y4(i)=r10t*ssi(i)*(r101r/zs(i)**2+r102rf/zs(i)**bsh5)/dd(i)
         psdep(i)=max(y4(i), 0.0E0)
          dd(i)=y3(i)*(rn20a*y3(i)-rn20b)+rn10c*tair(i)/esi(i)
          y2(i)=r191r/zg(i)**2+r192rf/zg(i)**bgh5
         pgdep(i)=R2ICE*max(r20t*ssi(i)*y2(i)/dd(i), 0.0e0)
c     ******************************************************************
          if (ice2 .eq. 0) then
            y1(i)=min(psdep(i)+pgdep(i), dd1(i))
           pgdep(i)=y1(i)-psdep(i)
          else
           psdep(i)=min(psdep(i), dd1(i))
           y1(i)=psdep(i)
          endif
         pt(i)=pt(i)+ascp*y1(i)            !LH
         rlh(i,k)=rlh(i,k)+ascp*y1(i)        !LH
         qv(i)=qv(i)-y1(i)
         qs(i)=qs(i)+psdep(i)
         qg(i)=qg(i)+pgdep(i)
        endif
c* 23 * ern : evaporation of qr                                   **23**
        ern(i)=0.0
        if (qr(i) .gt. 0.0) then
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
          pt(i)=pt(i)-avcp*ern(i)        !LH
          rlh(i,k)=rlh(i,k)-avcp*ern(i)    !LH
          qv(i)=qv(i)+ern(i)
          qr(i)=qr(i)-ern(i)
         endif
c* 30 * pmltg : evaporation of melting qg                         **30**
c* 33 * pmlts : evaporation of melting qs                         **33**
        pmlts(i)=0.0
        pmltg(i)=0.0
         tair(i)=(pt(i)+tb0)*pi0
        if (tair(i) .ge. t0) then
           rtair(i)=1./(t0-c358)
            y2(i)=exp( c172-c409*rtair(i) )
           esw(i)=c610*y2(i)
           qsw(i)=rp0*y2(i)
           ssw(i)=1.-(qv(i)+qb0)/qsw(i)
           dm(i)=qsw(i)-qv(i)-qb0
           rsub1(i)=cv409*qsw(i)*rtair(i)*rtair(i)
          dd1(i)=max(dm(i)/(1.+rsub1(i)),0.0E0)
            y3(i)=1./tair(i)
           dd(i)=y3(i)*(rn30a*y3(i)-rn10b)+rn10c*tair(i)/esw(i)
           y1(i)=r30t*ssw(i)*(r191r/zg(i)**2+r192rf/zg(i)**bgh5)/dd(i)
          pmltg(i)=R2ICE*min(qg(i),max(y1(i),0.0E0))
           y1(i)=r33t*ssw(i)*(r331r/zs(i)**2+r332rf/zs(i)**bsh5)/dd(i)
          pmlts(i)=min(qs(i),max(y1(i),0.0E0))
          if (ice2 .eq. 1) then
             y1(i)=min(pmlts(i), dd1(i))
            pmlts(i)=y1(i)
          else
             y1(i)=min(pmltg(i)+pmlts(i),dd1(i))
            pmltg(i)=y1(i)-pmlts(i)
          endif
          pt(i)=pt(i)-ascp*y1(i)       ! LH
          rlh(i,k)=rlh(i,k)-ascp*y1(i)   ! LH
          qv(i)=qv(i)+y1(i)
          qs(i)=qs(i)-pmlts(i)
          qg(i)=qg(i)-pmltg(i)
        endif
c        if (qv(i)+qb0 .le. 0.) qv(i)=-qb0
        if(qc(i).lt.cmin) qc(i)=0.0
        if(qr(i).lt.cmin) qr(i)=0.0
        if(qi(i).lt.cmin) qi(i)=0.0
        if(qs(i).lt.cmin) qs(i)=0.0
        if(qg(i).lt.cmin) qg(i)=0.0
       dpt(i,k)=pt(i)
       dqv(i,k)=qv(i)
       qcl(i,k)=qc(i)
       qrn(i,k)=qr(i)
       qci(i,k)=qi(i)
       qcs(i,k)=qs(i)
       qcg(i,k)=qg(i)

c       q1_hyd(i,k)=q1_hyd(i,k)+pt(i)-tttbud(i)
c       q2_hyd(i,k)=q2_hyd(i,k)+qv(i)-qqqbud(i)
c
c       if(isec.eq.isec/ibudsec*ibudsec) then
c          q1a_hyd(i,k) = q1_hyd(i,k) / rbud
c          q2a_hyd(i,k) = q2_hyd(i,k) / rbud
c          q1_hyd(i,k)=0.
c          q2_hyd(i,k)=0.
c       endif


  300  continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       scc=0.
       see=0.
           sddd=0.
           ssss=0.
           shhh=0.
           sccc=0.
           smmm=0.
           sfff=0.
       do 105 i=it1,it2
         dd(i)=max(-cnd(i)-cnd1(i), 0.e0)
         cnd(i)=max(cnd(i)+cnd1(i), 0.e0)
         dd1(i)=max(-dep(i)-dep1(i), 0.e0)
         dep(i)=max(dep(i)+dep1(i), 0.e0)
c
         sddd=sddd+dep(i)+pint(i)+psdep(i)+pgdep(i)
         ssss=ssss+dd1(i)
         shhh=shhh+rsw(i,k)*dt
         sccc=sccc+rlw(i,k)*dt
         smmm=smmm+psmlt(i)+pgmlt(i)+pimlt(i)
         sfff=sfff+(psacw(i)+piacr(i)+psfw(i)+pgfr(i)
     1        +dgacw(i)+dgacr(i)+psacr(i))*dt-qracs(i)+pihom(i)+pidw(i)
c
         scc=scc+cnd(i)
  105    see=see+dd(i)+ern(i)
        sc(k)=scc+sc(k)
        se(k)=see+se(k)
c
        s_dep(k)=s_dep(k)+sddd
        s_sub(k)=s_sub(k)+ssss
        s_qrs(k)=s_qrs(k)+shhh
        s_qrl(k)=s_qrl(k)+sccc
        s_mel(k)=s_mel(k)+smmm
        s_frz(k)=s_frz(k)+sfff
c
        ENDIF
c
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
          y1(i)=qc(i)+qr(i)+qi(i)+qs(i)+qg(i)
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
          a2=0.0
          a3=0.0
          a4=0.0
          a5=0.0
          a6=0.0
          a7=0.0
          a8=0.0
          a9=0.0
          a10=0.0
          a11=0.0
          a12=0.0
          a13=0.0
          a14=0.0
          a15=0.0
          a16=0.0
          a17=0.0
          a18=0.0
          a19=0.0
          a20=0.0
          a21=0.0
          a22=0.0
          a23=0.0
          a24=0.0
          a25=0.0
          a26=0.0
          a27=0.0
          a28=0.0
          a29=0.0
          a30=0.0
          a31=0.0
          a32=0.0
          a33=0.0
          a34=0.0
          a35=0.0
          a36=0.0
          a37=0.0
          a38=0.0
          a39=0.0
          a40=0.0
          a41=0.0
          a42=0.0
C
          TMODSUM = 0.
          QMODSUM = 0.
C

          do 30 i=it1,it2
          if(ibz(i,mt).eq.1) then
          if(kc.eq.4)then

c          bs11=avc*(cnd(i)-ern(i)-dd(i))*ndt_stat/d2t
          bs12=asc*(dep(i)-dd1(i)+psdep(i)+pgdep(i)+pint(i))
     1       *ndt_stat/d2t
c          bs13=afc*(psmlt(i)+pgmlt(i)+pimlt(i))*ndt_stat/d2t
c          bs14=afc*(psacw(i)+piacr(i)+psfw(i)+pgfr(i)+dgacw(i)
c     1       +dgacr(i)+psacr(i))*ndt_stat-afc*qracs(i)*ndt_stat/d2t
c     2       +afc*(pihom(i)+pidw(i))*ndt_stat/d2t
          bs15=asc*(pmltg(i)+pmlts(i))*ndt_stat/d2t

c          ceds1i(i,k,mt)=ceds1i(i,k,mt)+bs11+bs12+bs14-bs13-bs15
c          ceds2i(i,k,mt)=ceds2i(i,k,mt)-bs11/avc
     1                     -(bs12-bs15)/asc
          endif
          sww=sww+dm(i)
          scc=scc+cnd(i)
C
          TMODSUM = TMODSUM + TMODO(I,K)
          QMODSUM = QMODSUM + QMODO(I,K)
C
          see=see+dd(i)+ern(i)
c          a1=a1+pwacs(i)
          a1=a1+pihom(i)+pidw(i)
          a2=a2+pint(i)
          a3=a3+pgfr(i)
          a4=a4+psaut(i)
          a5=a5+psaci(i)
          a6=a6+psacw(i)
          a7=a7+praci(i)
          a8=a8+piacr(i)
          a9=a9+praut(i)
          a10=a10+pracw(i)
          a11=a11+psfw(i)
          a12=a12+psfi(i)
          a13=a13+pgacs(i)+dgacs(i)
          a14=a14+dgacw(i)
          a15=a15+dgaci(i)
          a16=a16+dgacr(i)
          a17=a17+pmltg(i)
          a18=a18+dep(i)
          a19=a19+pracs(i)
          a20=a20+psacr(i)
          a21=a21+pmlts(i)
          a22=a22+psmlt(i)
          a23=a23+pgmlt(i)
          a24=a24+psdep(i)
          A25=A25+PIMLT(I)
          a26=a26+pgdep(i)
          a27=a27+dd1(i)
          a28=a28+praci(i)*dlt3(i)
          a29=a29+piacr(i)*dlt3(i)
          a30=a30+psacr(i)*dlt2(i)
          a31=a31+qracs(i)*rdt
          a32=a32+psacw(i)*dlt4(i)
          a33=a33+qcl(i,k)
          a34=a34+qrn(i,k)
          a35=a35+qci(i,k)
          a36=a36+qcs(i,k)
          a37=a37+qcg(i,k)
          a38=a38+ern(i)
          a39=a39+wgacr(i)
          a40=a40+qsacw(i)
          a41=a41+dda(i)
          a42=a42+ddb(i)
          endif
   30   continue
        smf0(k,mt,kc)=sww+smf0(k,mt,kc)
        coc(k,mt,kc)=scc+coc(k,mt,kc)
        coe(k,mt,kc)=see+coe(k,mt,kc)
c
        t_bm(k,mt,kc) = TMODSUM + t_bm(k,mt,kc)
        q_bm(k,mt,kc) = QMODSUM + q_bm(k,mt,kc)
c
        thom(k,mt,kc)=thom(k,mt,kc)+a1
        tdw(k,mt,kc)=tdw(k,mt,kc)+a2
        tmlt(k,mt,kc)=tmlt(k,mt,kc)+a3
        saut(k,mt,kc)=saut(k,mt,kc)+a4
        saci(k,mt,kc)=saci(k,mt,kc)+a5
        sacw(k,mt,kc)=sacw(k,mt,kc)+a6
        raci(k,mt,kc)=raci(k,mt,kc)+a7
        tacr(k,mt,kc)=tacr(k,mt,kc)+a8
        raut(k,mt,kc)=raut(k,mt,kc)+a9
        racw(k,mt,kc)=racw(k,mt,kc)+a10
        sfw(k,mt,kc)=sfw(k,mt,kc)+a11
        sfi(k,mt,kc)=sfi(k,mt,kc)+a12
        gacs(k,mt,kc)=gacs(k,mt,kc)+a13
        gacw(k,mt,kc)=gacw(k,mt,kc)+a14
        gaci(k,mt,kc)=gaci(k,mt,kc)+a15
        gacr(k,mt,kc)=gacr(k,mt,kc)+a16
        gwet(k,mt,kc)=gwet(k,mt,kc)+a17
        gaut(k,mt,kc)=gaut(k,mt,kc)+a18
        racs(k,mt,kc)=racs(k,mt,kc)+a19
        sacr(k,mt,kc)=sacr(k,mt,kc)+a20
        gfr(k,mt,kc)=gfr(k,mt,kc)+a21
        smlt(k,mt,kc)=smlt(k,mt,kc)+a22
        gmlt(k,mt,kc)=gmlt(k,mt,kc)+a23
        sdep(k,mt,kc)=sdep(k,mt,kc)+a24
        ssub(k,mt,kc)=ssub(k,mt,kc)+a25
        gsub(k,mt,kc)=gsub(k,mt,kc)+a26
        pern(k,mt,kc)=pern(k,mt,kc)+a27
        d3ri(k,mt,kc)=d3ri(k,mt,kc)+a28
        d3ir(k,mt,kc)=d3ir(k,mt,kc)+a29
        d2sr(k,mt,kc)=d2sr(k,mt,kc)+a30
        d2rs(k,mt,kc)=d2rs(k,mt,kc)+a31
        gdry(k,mt,kc)=gdry(k,mt,kc)+a32
        erns(k,mt,kc)=erns(k,mt,kc)+a38
        wgrs(k,mt,kc)=wgrs(k,mt,kc)+a39
        qsws(k,mt,kc)=qsws(k,mt,kc)+a40
        srsw(k,mt,kc)=srsw(k,mt,kc)+a41
        srlw(k,mt,kc)=srlw(k,mt,kc)+a42
        qc0(k,mt,kc)=a33*ril2
        qr0(k,mt,kc)=a34*ril2
        qi0(k,mt,kc)=a35*ril2
        qs0(k,mt,kc)=a36*ril2
        qg0(k,mt,kc)=a37*ril2
        sqc0(k,mt,kc)=sqc0(k,mt,kc)+qc0(k,mt,kc)
        sqr0(k,mt,kc)=sqr0(k,mt,kc)+qr0(k,mt,kc)
        sqi0(k,mt,kc)=sqi0(k,mt,kc)+qi0(k,mt,kc)
        sqs0(k,mt,kc)=sqs0(k,mt,kc)+qs0(k,mt,kc)
        sqg0(k,mt,kc)=sqg0(k,mt,kc)+qg0(k,mt,kc)
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
C    DEPOSITION:    DEP(IJ)+PSDEP(IJ)+PGDEP(IJ)+PINT(IJ)
C    SUBLIMATION:   DD1(IJ)+PMLTS(IJ)+PMLTG(IJ)
C    MELTING:       PSMLT(IJ)+PGMLT(IJ)+PIMLT(IJ)+QRACS(IJ)
C    FREEZING:      PIHOM(IJ)+PIDW(IJ)+PSACW(IJ)+PSFW(IJ)+DGACW(IJ)
C                   +PIACR(IJ)+DGACR(IJ)+PSACR(IJ)+PGFR(IJ)
C    MASS FLUX:     DM(IJ)
C    CLOUD WATER:   QC(IJ)
C    RAIN:          
C    CLOUD ICE
C    SNOW
C    HAIL/GRAUPEL:
C----------------------------------------------------------------------
       DO 42 IJ=it1,it2
        CND(IJ)=CND(IJ)*RFT
        ERN(IJ)=(ERN(IJ)+DD(IJ))*RFT
        Y1(IJ)=(DEP(IJ)+PSDEP(IJ)+PGDEP(IJ)+PINT(IJ))*RFT
        Y2(IJ)=(DD1(IJ)+PMLTS(IJ)+PMLTG(IJ))*RFT
        Y3(IJ)=(PSMLT(IJ)+PGMLT(IJ)+PIMLT(IJ)+QRACS(IJ))*RFT
        Y4(IJ)=(PSACW(IJ)+PSFW(IJ)+DGACW(IJ)+
     1          PIACR(IJ)+DGACR(IJ)+PSACR(IJ)+PGFR(IJ)+(pihom(ij)+
     2          pidw(ij))*rdt)*RFT

        Y5(IJ)=DM(IJ)
        A1=1.E6*R00*QR(IJ)                                                  
        A2=1.E6*R00*QS(IJ)                                                  
        A3=1.E6*R00*QG(IJ)                                                  
        A11=UCOR*(MAX(1.E-5,A1))**1.75
        A22=UCOS*(MAX(1.E-5,A2))**1.75
        A33=UCOG*(MAX(1.E-5,A3))**1.75
        ZDRY=MAX(1.,A11+A22+A33)
        DBZ(IJ)=10.*LOG10(ZDRY)
        IF (TAIR(IJ).GE.273.16) THEN
          A44=A11+UWET*(A22+A33)**.95
          ZWET=MAX(1.,A44)
          DBZ(IJ)=10.*LOG10(ZWET)
        ENDIF
        DBZ(IJ)=DBZ(IJ)
        QC(IJ)=QC(IJ)
        QR(IJ)=QR(IJ)
        QI(IJ)=QI(IJ)
        QS(IJ)=QS(IJ)
        QG(IJ)=QG(IJ)

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
          S11(IWW,K)=S11(IWW,K)+Y1(IJ)
          S12(IWW,K)=S12(IWW,K)+Y2(IJ)
          S13(IWW,K)=S13(IWW,K)+Y3(IJ)
          S14(IWW,K)=S14(IWW,K)+Y4(IJ)
          S15(IWW,K)=S15(IWW,K)+Y5(IJ)
          S16(IWW,K)=S16(IWW,K)+QC(IJ)
          S17(IWW,K)=S17(IWW,K)+QR(IJ)
          S18(IWW,K)=S18(IWW,K)+QI(IJ)
          S19(IWW,K)=S19(IWW,K)+QS(IJ)
          S20(IWW,K)=S20(IWW,K)+QG(IJ)
          S21(IWW,K)=S21(IWW,K)+DBZ(IJ)
          SCNT(IWW,K)=SCNT(IWW,K)+1.
        ENDIF
        IF (ICS(IJ,2) .EQ. 1 .AND. JWW .GE. 1) THEN
          SN9(jww,K)=SN9(jww,K)+CND(IJ)
          SN10(jww,K)=SN10(jww,K)+ERN(IJ)
          SN11(jww,K)=SN11(jww,K)+Y1(IJ)
          SN12(jww,K)=SN12(jww,K)+Y2(IJ)
          SN13(jww,K)=SN13(jww,K)+Y3(IJ)
          SN14(jww,K)=SN14(jww,K)+Y4(IJ)
          SN15(jww,K)=SN15(jww,K)+Y5(IJ)
          SN16(jww,K)=SN16(jww,K)+QC(IJ)
          SN17(jww,K)=SN17(jww,K)+QR(IJ)
          SN18(jww,K)=SN18(jww,K)+QI(IJ)
          SN19(jww,K)=SN19(jww,K)+QS(IJ)
          SN20(jww,K)=SN20(jww,K)+QG(IJ)
          SN21(jww,K)=SN21(jww,K)+DBZ(IJ)
          SNCNT(jww,K)=SNCNT(jww,K)+1.
        ENDIF
        IF (ICS(IJ,3) .EQ. 1 .AND. IWW .GE. 1) THEN
          SS9(IWW,K)=SS9(IWW,K)+CND(IJ)
          SS10(IWW,K)=SS10(IWW,K)+ERN(IJ)
          SS11(IWW,K)=SS11(IWW,K)+Y1(IJ)
          SS12(IWW,K)=SS12(IWW,K)+Y2(IJ)
          SS13(IWW,K)=SS13(IWW,K)+Y3(IJ)
          SS14(IWW,K)=SS14(IWW,K)+Y4(IJ)
          SS15(IWW,K)=SS15(IWW,K)+Y5(IJ)
          SS16(IWW,K)=SS16(IWW,K)+QC(IJ)
          SS17(IWW,K)=SS17(IWW,K)+QR(IJ)
          SS18(IWW,K)=SS18(IWW,K)+QI(IJ)
          SS19(IWW,K)=SS19(IWW,K)+QS(IJ)
          SS20(IWW,K)=SS20(IWW,K)+QG(IJ)
          SS21(IWW,K)=SS21(IWW,K)+DBZ(IJ)
          SSCNT(IWW,K)=SSCNT(IWW,K)+1.
        ENDIF
        IF(ICS(IJ,3) .EQ. 1 .AND. jww.GE.1) THEN
          SSN9(jww,K)=SSN9(jww,K)+CND(IJ)
          SSN10(jww,K)=SSN10(jww,K)+ERN(IJ)
          SSN11(jww,K)=SSN11(jww,K)+Y1(IJ)
          SSN12(jww,K)=SSN12(jww,K)+Y2(IJ)
          SSN13(jww,K)=SSN13(jww,K)+Y3(IJ)
          SSN14(jww,K)=SSN14(jww,K)+Y4(IJ)
          SSN15(jww,K)=SSN15(jww,K)+Y5(IJ)
          SSN16(jww,K)=SSN16(jww,K)+QC(IJ)
          SSN17(jww,K)=SSN17(jww,K)+QR(IJ)
          SSN18(jww,K)=SSN18(jww,K)+QI(IJ)
          SSN19(jww,K)=SSN19(jww,K)+QS(IJ)
          SSN20(jww,K)=SSN20(jww,K)+QG(IJ)
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
       sq(k)=sq(k)+qcl(i,k)+qrn(i,k)+qci(i,k)+qcg(i,k)+qcs(i,k)
       sqc(k)=sqc(k)+qcl(i,k)
       sqr(k)=sqr(k)+qrn(i,k)
       sqi(k)=sqi(k)+qci(i,k)
       sqs(k)=sqs(k)+qcs(i,k)
       sqg(k)=sqg(k)+qcg(i,k)
  465 continue
      do 470 k=2,kles
       y1(k)=y1(k)*ril2
       y2(k)=y2(k)*ril2
       sq(k)=sq(k)*ril2
       sqc(k)=sqc(k)*ril2
       sqr(k)=sqr(k)*ril2
       sqi(k)=sqi(k)*ril2
       sqs(k)=sqs(k)*ril2
       sqg(k)=sqg(k)*ril2
       stqc(k)=stqc(k)+sqc(k)
       stqr(k)=stqr(k)+sqr(k)
       stqi(k)=stqi(k)+sqi(k)
       stqs(k)=stqs(k)+sqs(k)
       stqg(k)=stqg(k)+sqg(k)
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
         QCI1(I,K)=QCI1(I,K)+EPS*QCI(I,K)
         QCS1(I,K)=QCS1(I,K)+EPS*QCS(I,K)
         QCG1(I,K)=QCG1(I,K)+EPS*QCG(I,K)
c       ELSE
c          DPT1(I,K)=DPT(I,K)
c          DQV1(I,K)=DQV(I,K)
c          QCL1(I,K)=QCL(I,K)
c          QRN1(I,K)=QRN(I,K)
c          QCI1(I,K)=QCI(I,K)
c          QCS1(I,K)=QCS(I,K)
c          QCG1(I,K)=QCG(I,K)
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
      call boundy (qci,qci1)
      call boundy (qcs,qcs1)
      call boundy (qcg,qcg1)
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
        FD(K)=(TA(K)-TA(K-1))*RHO1(K)*.5
        FE(K)=(QA(K)-QA(K-1))*RHO1(K)*.5
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
