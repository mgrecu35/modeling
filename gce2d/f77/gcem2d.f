C              if using cray, compile with -eZ flag
C              uncomment following line if using sgi, add comment if using cray!

      PROGRAM COARE_GATE

C  SC0 changed to 1411 for TOGA WKSHOP #2
C
C 1/16/97 F(x,t)->  ths - sfc theta; qs - sfc sat mr; ts - sfc temp.
C   F(t)-> thsdt - 6hr sfc. theta tend.; qsdt -6 hr sfc sat m.r. tend.
C                         tsdt - 6 hr sfc. temp. tend.
C   F(x)-> tairsfc - sfc. air temp; qairsfc - sfc. sat m.r.; 
C          thairsfc - sfc. theta
C
c            July 27 2001 (Tao)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    *****   THIS CODE HAS BETTS AND MILLER SCHEME                     C
C    *****   THIS CODE CAN RUN ONLY 2ND OR 4TH ORDER ADVECTION SCHEMES C
C    *****   FILENAME=TAO2D                                            C
C    *****   CYCLIC LATERAL BOUNDARY CONDITION                         C
C    *****   FOUR ADDITIONAL OPTIONS                                   C
C    *****   IMLIFTING = 1  FOLLOWED SOONG & OGURA                     C
C    *****   IMLIFTING = 0  TRADITIONAL CYCLIC B. C. MODEL             C
C    *****   LIPPS = 1      FOLLOWED LIPPS & HEMLER'S MODIFICATION     C
C    *****   LIPPS = 0      ORIGINAL ANELASTIC SYSTEM                  C
C    *****   IJKADV = 1  POSITIVE DEFINITE ADVECTION SCHEME            C
C    *****   IJKADV = 0  2ND AND/OR 4TH ORDER ADVECTION SCHEME         C
C    *****   ISNGTAO = 1  D(THEA)/DZ PREDICTED AS SOONG AND TAO        C
C    *****   ISNGTAO = 0  D(THEA)/DZ PRECRIBED                         C
C    *****   IWBAR = 1    IWBAR IS CONSTANT AS SUI ET AL               C
C    *****   IWBAR = 0    IWBAR VARIES AS IMPOSED  TAO AND SOONG       C
C----------------------------------------------------------------------C
C    *****                                                             C
C    *****   YOU HAVE TO HAVE FOLLOWING FILES TO RUN THIS PROGRAM:     C
C    *****       DATAFILE --- PARAMETER CONTROL FILE                   C
C    *****       cah.dat, cai.dat,                                     C
C    *****       co2.tran3, h2o.tran3, o3.tran3                        C
C    *****       trp.dat  --- INCLUDES FILES FOR THE RADIATION PART    C
C    *****       coare_st_data4_225 --- Sui's INITIAL SOUNDING DATA    C
C    *****                              used in base routine           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c                                                                      c
c                                                                      c
C    NEED TO CHANGE NXX WHEN NX IS CHANGED                             C
c                                                                      c
c                                                                      c
c     ****   slab-symmetric cloud ensemble model (starts 3/1/76)       c
c     ****   second order flux form (starts 3/9/78)                    c
c     ****   bob"s diffusion coefficient (starts 3/9/78)               c
c     ****   stretched in z-direction (starts 10/8/80)                 c
c     ****   for cyber 205 vectorized processor (starts 12/18/82)      c
c     ****   open lateral boundary condition  (starts 3/18/84)         c
c     ****   3-classes ice-phase  (starts 12/18/85)                    c
c     ****   fourth order flux form (starts 6/9/86)                    c
c     ****   cray (starts 9/9/90)                                      c
c     ****   changes made for a long integration (12/05/91)            c
c              triple integratin time to 18 hours for each run         c
c              imp,imx,imx1 are tripled, iwp increased to 3600         c
C     ****   portable random numbers generator from "Num.Rec."(1986)   c
C              8/Mar/1994                                              c 
C     ****   change variables AMIN to AMINUT as needed                 c
c                                                                      c
c     modify I/O subroutines ----- April 28, 94                        c
c1    subroutine wtap(x,itape) divide into wtap(x,itape) &             c
c1       wtap_sngl(x,itape)                                            c
c1    subroutine wtapt(itape)  divide into  wtapt(itape) &             c
c1       wtapt_sngl(itape)                                             c 
c2    subroutine wstart (itape) unchanged                              c
c3    the following subroutines are all modified to output             c
c3       real*4 variables                                              c
c3          wtapri, wtaps(x,n), wtap3(x), wtapp                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccshie
c     parameter (NX 66,NZ=43,NT2304,ITT=244,lay=88,NXI=NX-2,nnt=481) ! 8d rain
c     parameter (NX 66,NZ=43,NT=2880,ITT=244,lay=88,NXI=NX-2,nnt=481) !20d rain
c     parameter (NX 258,NZ=43,NT=2880,ITT=244,lay=88,NXI=NX-2,nnt=481) ! 4/18/01
c     parameter (NX 130,NZ=43,NT=2880,ITT=244,lay=88,NXI=NX-2,nnt=481) ! 4/18/01
c     parameter (NX 514,NZ=43,NT=2880,ITT=244,lay=88,NXI=NX-2,nnt=481) ! 4/18/01
c     parameter (NX 258,NZ=43,NT 5760,ITT=244,lay=88,NXI=NX-2,nnt=481) !6/7/01
      parameter (NX=514,NZ=43,NT=2880,ITT=244,lay=88,NXI=NX-2,nnt=481) !6/7/01
      parameter (nx7=7*nx,nt3=3*nt,itt3=3*itt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      common/radflux/rflux(nxi/2,7)
      COMMON/SFLUXS/ SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      COMMON/PBLTEN/ UTEN(NX,NZ),VTEN(NX,NZ),TTEN(NX,NZ),
     1      QTEN(NX,NZ),QCTEN(NX,NZ),QCITEN(NX,NZ),PBLDT,ATIME
      common/itoga/ itoga,ISFC,ICE,ICE2
      common/dinrad/p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),thairsf(nx)
     1,     pairsfc(nx)
      common/sfcri/ ri180(nx),riold180(nx),iricont
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/mbudget/ acoc(nz),acoe(nz),agaut(nz),asdep(nz),agsub(nz),
     1   apern(nz),aqls(nz),altrans(nz),actrans(nz),total_mq(nz),
     2   agwet(nz),afgr(nz),avelw(nz),avesw(nz),amelt(nz),aother(nz),
     3   total_mt(nz)
c     time sequence of model domain averaged t and q
      common/tqave/ tavet(nz),qavet(nz)
      common/tbudget/ avett(nnt),avetq(nnt),avet950(nnt)
      common/bsfc/ tsfc_1(nx),qsfc_1(nx)
      common/gbs/ tlsw(nz),qlsw(nz),ttlsw(nz),tqlsw(nz),smfff(nz),
     1   smffff(nz),tcof(nz),textra(nz),qextra(nz)
      common/gbs11/ tlsw1(nz),qlsw1(nz),ttlsw1(nz),tqlsw1(nz) 
      common/timestat/ ndt_stat,itime_ave,mmave
      common/tvertical/ denave
      common/tbudget1/ t_sfcq(nnt),t_sfct(nnt),t_ub(nnt),t_microq(nnt),
     1   t_largeq(nnt),t_larget(nnt),t_microt(nnt)
c      dimension total_mqt(nz),total_mtt(nz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/msui/ isounding,isui,ifine,idelz
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      COMMON/IPTIONR/ IRADAVE
      COMMON/IPTIONR1/ IOPCLOUD
      common/iice/ new_ice_sat
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
      common/bstart/ n,isec,nran
      common/rbstart/ sec,aminut,rdt
      common/bcor/ irf,iadvh,irfg,idq,ismg
      common/bcorr/ id
      common/iceopt/ ice913,ilif
      common/ilin/ lin
      common/bcor1/ fcor
      common/bsize/ it1,it2,kt1,kt2
      common/bnumer/ kref
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
      common/b1u/ u(nx,nz)
      common/b1v/ v(nx,nz)
      common/b1w/ w(nx,nz)
      common/b1a/ ak(nx,nz)
      common/b2u/ uu1(nx,nz)
      common/b2v/ vv1(nx,nz)
      common/b2w/ ww1(nx,nz)
      common/b2a/ ak1(nx,nz)
      common/bsat/ xxw(nx,nz)
      common/bsat1/ aaa(nx,nz)
      common/badv/ dfc(nx,nz)
      COMMON/BTV/ VTP(nx,nz)
      common/dumuw/ umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      COMMON/BW/ XXX(NX,NZ),DF0(NX,NZ),TRAH(NX,NZ),TRAV(NX,NZ)
      common/slwave/rsw(nx,nz),rlw(nx,nz)
      common/bw1/ rst(nt,nxi),pcltop(nt,nxi),pclbot(nt,nxi)
      common/bw2/ rsv(nt,nxi),RSVAR(nt,nxi)
      common/bw3/ rwmax(nt3),rwmin(nt3),anv(nt3),cnv(nt3),snp(nt3)

      common/b4/ tbsk(nz),bskt(nz),bskt2(nz),bskm(nz),bskm4(nz),
     $   bsit(nz),bsit2(nz),bsim(nz),bsim4(nz)
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b8/ smf(nz),smu(nz),smd(nz),smn(nz),stf(nz),stu(nz),
     $   std(nz),stn(nz),sqf(nz),squ(nz),sqd(nz),sqn(nz),sdt(nz),
     $   sdq(nz),stv(nz),
     $   stf1(nz),stu1(nz),std1(nz),stn1(nz),sqf1(nz),squ1(nz),
     $   sqd1(nz),sqn1(nz),stf2(nz),stu2(nz),std2(nz),stn2(nz),sqf2(nz),
     $   squ2(nz),sqd2(nz),sqn2(nz),cld(nz),cldu(nz),cldd(nz),cldn(nz)
      common/b9/ sut1(nz),suc1(nz),sun1(nz),suu1(nz),sud1(nz),aub(nz),
     1   ssu(nz),sut2(nz),suc2(nz),sun2(nz),suu2(nz),sud2(nz),auc1(nz),
     2   aun1(nz),auu1(nz),aud1(nz),auc2(nz),aun2(nz),auu2(nz),aud2(nz),
     3   svt1(nz),svc1(nz),svn1(nz),svu1(nz),svd1(nz),avb(nz),
     4   ssv(nz),svt2(nz),svc2(nz),svn2(nz),svu2(nz),svd2(nz),avc1(nz),
     5   avn1(nz),avu1(nz),avd1(nz),avc2(nz),avn2(nz),avu2(nz),avd2(nz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b66b/ s_dep(nz),s_sub(nz),s_qrs(nz),s_qrl(nz),s_mel(nz),
     1   s_frz(nz)
      common/bb6/ tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/sfcten/ tsdt,qsdt,thsdt,psdt
      COMMON/Q1Q2Z/ Q1Z_6H(NZ,ITT),Q1ZT(NZ),Q2Z_6H(NZ,ITT),Q2ZT(NZ)
      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/bls3/ factor_nuding
      common/bls4/ ubi(nz),vbi(nz),ub_2h(nz,itt3),vb_2h(nz,itt3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/ba/ y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx7),
     1        tair(nx),qss(nx),qv(nx)
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
     $ tb0(nz,4),qb0(nz,4)
c      common/bstsi/ ceds1i(nx,nz,4),ceds2i(nx,nz,4),
c     1   tstfi(nx,nz,4),tsqfi(nx,nz,4),rsli(nx,nz,4)
      common/bsts1/ tut1(nz,4,7),tut2(nz,4,7),tvt1(nz,4,7),tvt2(nz,4,7),
     $ tstf(nz,4,7),tstf1(nz,4,7),tstf2(nz,4,7),tsqf(nz,4,7),
     $ tsqf1(nz,4,7),tsqf2(nz,4,7),tsqq(nz,4,7),tsqq1(nz,4,7),
     $ fcld(nz,4,7)
      COMMON/BSTS20/ OTHERT_ADD(NZ,4,7),OTHERQ_ADD(NZ,4,7)

      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNQH(NZ,4,7),
     1 SNQV(NZ,4,7),SNQD(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7)


      common/stls/ srsw(nz,4,7),srlw(nz,4,7)
      COMMON/BTT/ S1(16,NZ,14),SN1(5,NZ,14)
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
      common/brh1/ srro(nz),qrro(nz),sqc(nz),sqr(nz),sqi(nz),sqs(nz),
     1   sqg(nz),stqc(nz),stqr(nz),stqi(nz),stqs(nz),stqg(nz),tqc(nt),
     2   tqr(nt),tqi(nt),tqs(nt),tqg(nt)
      common/bch/ it(nx),iv(nt),ics(nx,4),ibz(nx,4)
      common/bch1/ rby(7)
      common/srflx/ sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)
      common/qcbf/ qcl_buffer(nx,nz)
      common/surface/ sun_4(nx,4)
      common/bls/y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/zlevel/cc1,cc2,z1(nz),z2(nz)
      common/o4/ a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
c
      common/cpsbm/ icps_bm,iexplicit
      common/sue/ ppress(nx,nz)
      common/sue1/ tmodo(nx,nz),qmodo(nx,nz),rainco(nx,1)
      common/sue3/ t_bm(nz,4,7),q_bm(nz,4,7)
      common/sue4/ ar_c(nx)
      common/add_f/ i_four
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    fv(nz),fv1(nz),qaa1(nz),RI120AR(NX)
      integer jv(nt),ised(12)
      character*1 aaa1(nx,nz),blk,dot,com,exm,pice,ast,ago,abt
c      real atsqf1(nz),atsqf2(nz),avefluxq(nz)
      real pl(lay),pa(lay)
      real bnxnz(nx,nz),bnz(nz),bnx4(nx,4),bntnxi(nt,nxi)
      real bnt(nt),bnt3(nt3),bnz47(nz,4,7),bnz4(nz,4),bnx7(nxi/2,7)
      real b16nz14(16,nz,14),b5nz14(5,nz,14),b16nz(16,nz),b5nz(5,nz)
c      real bnnt(nnt)
      real bnx(nx)
c      real bnxnz4(nx,nz,4)
C
CC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFF_VX(NX,NZ) :  V-TURBULENT IN X-COMPNT                        C
C     DIFF_VZ(NX,NZ) :  V-TURBULENT IN Z-COMPNT                        C
C     GRID_VX(NX,NZ) :  V-GRID SCALE TRANSPT IN X-COMPNT               C
C     GRID_VZ(NX,NZ) :  V-GRID SCALE TRANSPT IN Z-COMPNT               C
C     DIFF_NV(NX,NZ) :  V-NUMERICAL FILTER                             C
C     V_LARGE(NX,NZ) :  LARGE-SCALE FORCING IN V WIND                  C
C     PRE_V(NX,NZ)   :  HORIZONTAL PREES FORCING                       C
C     DT_VWIND(NX,NZ):  LOCAL TIME CHANGE TERM                         C
C     VB_MEAN(NZ,4)  :  MEAN V WIND                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFF_UX(NX,NZ) :  U-TURBULENT IN X-COMPNT                        C
C     DIFF_UZ(NX,NZ) :  U-TURBULENT IN Z-COMPNT                        C
C     GRID_UX(NX,NZ) :  U-GRID SCALE TRANSPT IN X-COMPNT               C
C     GRID_UZ(NX,NZ) :  U-GRID SCALE TRANSPT IN Z-COMPNT               C
C     DIFF_NU(NX,NZ) :  U-NUMERICAL FILTER                             C
C     U_LARGE(NX,NZ) :  LARGE-SCALE FORCING IN U WIND                  C
C     PRE_U(NX,NZ)   :  HORIZONTAL PREES FORCING                       C
C     DT_UWIND(NX,NZ):  LOCAL TIME CHANGE TERM                         C
C     UB_MEAN(NZ,4)  :  MEAN U WIND                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real            DIFF_VX(NX,NZ),DIFF_VZ(NX,NZ),GRID_VX(NX,NZ),
     1                GRID_VZ(NX,NZ),DIFF_NV(NX,NZ),V_LARGE(NX,NZ),
     2                PRE_V(NX,NZ),DT_VWIND(NX,NZ)
      COMMON/DEBUG_V/ DIFF_VX,DIFF_VZ,GRID_VX,
     1                GRID_VZ,DIFF_NV,V_LARGE,
     2                PRE_V,DT_VWIND

      real            DIFF_UX(NX,NZ),DIFF_UZ(NX,NZ),GRID_UX(NX,NZ),
     1                GRID_UZ(NX,NZ),DIFF_NU(NX,NZ),U_LARGE(NX,NZ),
     2                PRE_U(NX,NZ),DT_UWIND(NX,NZ)
      COMMON/DEBUG_U/ DIFF_UX,DIFF_UZ,GRID_UX,
     1                GRID_UZ,DIFF_NU,U_LARGE,
     2                PRE_U,DT_UWIND

      real             SDIFF_VX(NZ,4),SDIFF_VZ(NZ,4),SGRID_VX(NZ,4),
     1                 SGRID_VZ(NZ,4),SDIFF_NV(NZ,4),SV_LARGE(NZ,4),
     2                 SPRE_V(NZ,4),SDT_V(NZ,4),TN_VWIND(NZ,4)
      COMMON/DEBUG_SV/ SDIFF_VX,SDIFF_VZ,SGRID_VX,
     1                 SGRID_VZ,SDIFF_NV,SV_LARGE,
     2                 SPRE_V,SDT_V,TN_VWIND

      real             SDIFF_UX(NZ,4),SDIFF_UZ(NZ,4),SGRID_UX(NZ,4),
     1                 SGRID_UZ(NZ,4),SDIFF_NU(NZ,4),SU_LARGE(NZ,4),
     2                 SPRE_U(NZ,4),SDT_U(NZ,4),TN_UWIND(NZ,4)
      COMMON/DEBUG_SU/ SDIFF_UX,SDIFF_UZ,SGRID_UX,
     1                 SGRID_UZ,SDIFF_NU,SU_LARGE,
     2                 SPRE_U,SDT_U,TN_UWIND

      real        UB_MEAN(NZ,4),VB_MEAN(NZ,4)
      COMMON/SUV/ UB_MEAN,VB_MEAN


      real             SSDIFF_VX(NZ,4),SSDIFF_VZ(NZ,4),SSGRID_VX(NZ,4),
     1                 SSGRID_VZ(NZ,4),SSDIFF_NV(NZ,4),SSV_LARGE(NZ,4),
     2                 SSPRE_V(NZ,4),SSDT_V(NZ,4),STN_VWIND(NZ,4)
      COMMON/DEBUG_SVS/ SSDIFF_VX,SSDIFF_VZ,SSGRID_VX,
     1                  SSGRID_VZ,SSDIFF_NV,SSV_LARGE,
     2                  SSPRE_V,SSDT_V,STN_VWIND

      real             SSDIFF_UX(NZ,4),SSDIFF_UZ(NZ,4),SSGRID_UX(NZ,4),
     1                 SSGRID_UZ(NZ,4),SSDIFF_NU(NZ,4),SSU_LARGE(NZ,4),
     2                 SSPRE_U(NZ,4),SSDT_U(NZ,4),STN_UWIND(NZ,4)
      COMMON/DEBUG_SUS/ SSDIFF_UX,SSDIFF_UZ,SSGRID_UX,
     1                  SSGRID_UZ,SSDIFF_NU,SSU_LARGE,
     2                  SSPRE_U,SSDT_U,STN_UWIND
C
CC
C
      real cosz,hrl,rlat  ! 8/21/01 shie
      integer month,iday  ! 8/21/01 shie

      data rby/2.,1.,0.,0.,0.,-.5,-1./
      data ised/54368,54368,33214,43598,21368,11246,78246,62359,31246,
     1   29548,83124,98272/
      data blk,dot,com,exm,pice,ast,ago,abt/' ','.',',','R','I','*'
     1                                      ,'G','O'/
      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ICASE=1:  Sui et al. 94   iwbar=1                          
C     ICASE=0:  TOGA/GATE       iwbar=0                          
C     ICE = 0       WARM RAIN ONLY
C     ICE2=1        QI AND QS ONLY
C     IFINE=0 SUI ET AL (1994)  COARSE VERTICAL RESOLUTION             
C     IFINE=1 GRABOWSKI ET AL (1996) FINE VERTICAL RESOLUTION          
C     IJKADV=0       HIGHER ORDER (2ND OR 4TH ORDER) ADVECTION    
C     IJKADV=1       LOWERER ORDERBUT POSITIVE DEFINITE ADVECTION 
C     ILIF=0  (MODEL INTERCOMPARISON WORKSHOP)
C     IMLIFTING = 1  FOLLOWED SOONG & OGURA		       
C     IMLIFTING = 0  TRADITIONAL CYCLIC B. C. MODEL	       
C     INCAR = 0 for Goddard data set - different unit than NCAR
C     IRADAVE=1:  Include cloud effects but apply Qr and sfc fluxes uniformly over domain 
C     ISNGTAO=0      WBAR * D (THEA)/DZ IS CONSTANT OR SPECIFIED 
C     ISNGTAO=1      WBAR IS THE MAIN FORCING & LAPSE RATE VARIES SEE SOONG AND TAO 1980  
C     ITOGA=2       AERODYNAMIC FORMULA
C     ITOGA=0:      BLACKADAR PBL
C     ITOGA=1:      TOGACOARE SURFACE FLUX ROURINE
C     IUVBAR = 1     U/VBAR ARE AS IMPOSED-SPECIFIED SEE SOONG AND TAO 1980  
C     IUVBAR = 0     U/VBAR VARY (NO CONTROL) SEE TAO AND SOONG 1986 
C     IWBAR = 1      WBAR IS CONSTANT WITH TIME - SEE SUI ET AL 1994 
C     IWBAR = 0      WBAR VARIES AS IMPOSED TAO AND SOONG  SEE SOONG AND TAO 1980  
C     ISOUNDING = 1 USE SUI ET AL'S SOUNDING                           
C     ISOUNDING = 0 USE GRABOWSKI ET AL'S SOUNDING                     
C     ISUI=1  SUI ET AL (1994)  WBAR * [D(THEA)/DZ]T                   
C     ISUI=0  GRABOWSKI ET AL (1996) [WBAR * D(THEA)/DZ]CONST          
C     LIPPS = 1      FOLLOWED LIPPS & HEMLER'S MODIFICATION      
C     LIPPS = 0      ORIGINAL ANELASTIC SYSTEM		      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     it starts -6 h and all temporal forcing are zero -> t,q,u,v and w
C                    are constant 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     RUN 1: TOGA COARE					               C
C            IMLIFTING=1  ISMGTAO=1                                    C
C     RUN 2: TOGA COARE						       C
C            IMLIFTING=0  ISMGTAO=0                                    C
C     RUN 3: GATE						       C
C            IMLIFTING=1  ISMGTAO=1                                    C
C     RUN 4: GATE						       C
C            IMLIFTING=0  ISMGTAO=0                                    C
C     ALL 4 RUNS:  IWBAR=0 IUVBAR=1 LIPPS=1 IDOM=1 IJKADV=0            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c       ibudsec=900
c       rbud = 1.

       call timing(time1)
       OPEN( 3,FILE='restart',FORM='UNFORMATTED')
       OPEN( 4,FILE='data',FORM='UNFORMATTED')
       OPEN( 7,FILE='start',FORM='UNFORMATTED',STATUS='UNKNOWN' )
c
cccccccccccccccccccc read input parameters ccccccccccccccccccccccccccccc
cccshie1
      irs=0  ! 1st time to run
c     irs=1  ! restart runs
      call rinit (irs)
c     inialize all parameters used in uand v wind statistics
c
      CALL RINIT_UV
c
      icps_bm=0
      iexplicit=1
        if (icps_bm .eq. 0) iexplicit=1
      ICE2=0
      if (icps_bm .eq. 1) ICE2=1
C     call cps index, initiated with 0
      KCPS=0
cccshie by tao 5/3/01
      new_ice_sat = 1 !  call tao et al (1989) saturation technique twice
c     new_ice_sat=0
c
      fcor=0.
      irf=1
      irfg=0
      isfc=1
cccshie
      itoga=1  ! toga flux
c     itoga=2  ! aerodynamic
      ice=1
      icase=0
      ilif=0
      ijkadv=1

c     i_four=1: 4th order filter in z-direction; only for ijkadv=0
c     i_four=0: 2nd order filter in z-direction; only for ijkadv=0 
      i_four=0

c     lipps=1
      lipps=0 ! cccshie by tao 5/3/01
      idom=0
      iradave=0
c     iradave=1  ! 10/4/01 sensitivity test for paper
      isngtao=0
c     imlifting=1
      imlifting=0 ! cccshie by tao 5/3/01
      iwbar=0
c
ctao 12-3-98
      iuvbar=1
ctao 12-3-98
c
      ifine=0
      ISOUNDING=0
      isui=0
      incar=1
      iadvh1=4
      irad=1
      rlat=0.0
c local time at sesmex (around 117E,22.5N), e.g., 0Z(Universal Time Coordinate)
c                                                 --> 8am local time
c 0Z, May 6, 1998
c     hrl=8.0 ! local time ! 4/20,4/18,4/8/01
c     hrl=8.0+6. ! local time ! 4/20,4/18,4/8/01
c     hrl=8.0+12. ! local time, 4/4/01, 5/3/01, may14.12z-may24.12z, NX 258
c     hrl=8.0+18. ! local time, 4/5/01,  4/6/01 not that hrl > 24 ??
c     hrl=8.0+8.  ! local time, 4/6/01, 4/9/01, 4/18/01 move forward 8 hrs
c     hrl=8.0+14.  ! local time, 4/18/01 move forward 14 hrs
c     iday=6   ! 4/20 re-run May 6, for 1st 10 days, but dt=6
c     iday=16  ! 4/3/01 start from May 16, instead of May 6.
c     iday=18  ! 4/6/01 4/9/01 start from May 18, instead of May 6.
c     iday=18  ! 8/1/01 start from May 18, for 5/18-5/20
c     iday=17  ! 6/7/01 6/8/01 tao proposed, May 17-22.
c     iday=20  ! 6/29/01, 2nd 4-5 day run 5/17-5/22
c     iday=33  ! 6/7/01 tao proposed, June 2-11. (31+2=33)
c     iday=14  ! 5/3/01, may14.12z-may24.12z, NX 258
c     iday=20  ! 4/6/01 4/9/01 start from May 20, 4/18/01
c     iday=21  ! 4/18, 4/8/01 start from May 21, try to pass the troble time

c     hrl=8.0+12. ! local time, 8/21/01 062612-062900, NX 258
      hrl=8.0  ! local time ! 8/29/01 051800-052600, NX 514
c     iday=57  ! 8/21/01 shallow convec, 062612-062900 (31+26=57)
c     iday=44  ! 8/21/01 shallow convec, 061312-061412 (31+13=44)
c     iday=14  ! 8/21/01 shallow convec, 051412-051512 ( 0+14=14)
c     iday=27  ! 8/21/01 shallow convec, 052712-052812 ( 0+27=27)
      iday=18  ! 8/29/01 long term integ, 051800-052600 ( 0+18=18)
c     iday=33  ! 8/29/01 long term integ, 060200-061100 (31+ 2=33)

      month=5
      nrun=2
cccshie2
c     time=1.  ! dan's test on 3/30/01
c     time=3.*24.*60.  ! in minutes, i.e., 3 days
c     time=10.*24.*60.  ! in minutes
c     time=5.*24.*60.  ! in minutes, 4/8/01, try may21-may26
c     time=5.*24.*60.  ! in minutes, 4/9/01, try may18-may23
c     time=3.*24.*60.  ! in minutes, 6/28/01, try 1st 3 days for may17-may22
c     time=2.*24.*60.  ! in minutes, 8/1/01, try 1st 2 days for may18-may20
      time=8.*24.*60.  ! in minutes, 8/29/01, try 8 days for may18-may26
c     time=9.*24.*60.  ! in minutes, 8/29/01, try 9 days for june2-june11
c     time=2.5*24.*60.  ! in minutes, 8/21/01, try 2.5 days for 062612-062900
c     time=24.*60.  ! in minutes, 8/21/01, try 1.0 day for 061312-061412
c     time=12.*60.  ! in minutes, 8/21/01, try 1st .5  day for 061312-061412
c     time=6.*60.  ! in minutes, 8/1/01, try 1st 2 days for may18-may20
c     time=12.*60.  ! in minutes, 8/1/01, try 1st 2 days for may18-may20
c     time=5.*24.*60.  ! in minutes, 6/29/01, try 2nd 4-5 days for may17-may22
c     time=9.*24.*60.  ! in minutes, 6/7/01, try june2-june11
c     time=10.*24.*60.  ! in minutes, 4/18/01, try may18-may28, NX 258
c     time=10.*24.*60.  ! in minutes, 5/3/01, may14.12z-may24.12z, NX 258
c      time=20.*24.*60.  ! in minutes  ! run 10-20 days with irs=1
	print *,'integrate to time (min,hr,day) = ',time,time/60.,time/1440.
c     dt=12.0  ! 4/6/01 again, 8/21/01 for nx 258
c     dt=10.0  ! 5/4/01 by tao for NX 514
      dt=7.5   ! dx=1000 m,  8/29/01, try 8 days for may18-may26
c     dt=6.0   ! 4/5/01, 4/18 for nx 258, 5/3/01 , may14.12z-may24.12z, NX 258
      idelz=0
      istaton=1
      imp=3600

cccshie 8/29/01 double iwp to reduce output data file size that integration 
c  days can be longer than 8, 9 days (with NX=514).
c  In machine "della1", maximum data file size = 2,147,483,647 bytes.
c     iwp=900
      iwp=1800  ! cccshie 8/29/01 to reduce file size that job can be run

      imx=300
      imx1=1800
      iwtp=0
      istrt=21600
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (isfc .eq. 0) itoga=10
      if (irad .eq. 0) iradd=0
      fcorg=fcor/980.
      PBLDT=180.          !compute pbl process every IPBLDT (IN SEC)
      IPBLDT=PBLDT
      idpbl=0
      ircl=0
CTAO
      ice913=0
      iopcloud=1
c     iopcloud=1: snow is treated as cloud ice
      lin=0
      icerh=1
      icelin=0
        if (icerh .eq. 1) icelin=0
        if (icelin .eq. 1) icerh=0
CTAO
c      isfcave = 1:  Use mean U and V for flux calculation
       isfcave=0
       fact_qvs=1.
CTAO
c     ******   set parameters for explicit vectorization
      ised1=ised(nrun)
      do i=1,nx
       ibz(i,1)=0
       ibz(i,2)=0
       ibz(i,3)=0
       ibz(i,4)=0
      end do
      do i=1,nx
       iv(i)=0
       jv(i)=0
       it(i)=0
       ics(i,1)=1
       ics(i,2)=0
       ics(i,3)=0
       ics(i,4)=0
       y0(i)=1.
       suw(i)=0.
       svw(i)=0.
       swt(i)=0.
       swq(i)=0.
       ri180(i)=0.
       RI120AR(I)=0.
       ar_c(i)=0.
      end do


      DO K=1,NZ
         Q1ZT(K)=0.
         Q2ZT(K)=0.
      DO I=1,ITT
         Q1Z_6H(K,I)=0.0
         Q2Z_6H(K,I)=0.0
      ENDDO
      ENDDO

      do k=1,nz
       smffff(k)=0.
       smfff(k)=0.
       tlsw(k)=0.
       qlsw(k)=0.
       ttlsw(k)=0.
       tqlsw(k)=0.
       tlsw1(k)=0.
       qlsw1(k)=0.
       ttlsw1(k)=0.
       tqlsw1(k)=0.
       tcof(k)=0.
c
       QAA1(K)=0.
       fv(k)=1.
       fv1(k)=1.
       tls(k)=0.0
       qls(k)=0.0
       UBT(K)=0.
       VBT(K)=0.
       WBT(K)=0.
       Q1T(K)=0.
       Q2T(K)=0.
       Q1ZT(K)=0.
       Q2ZT(K)=0.
      end do
       tsdt=0.
       qsdt=0.
       thsdt=0.
       psdt=0.
       aco=0.
       aan=0.
       aco1=0.
       aan1=0.
ccccc  ck and ce are parameters for sub-grid-scale turbulence
      ck=.2
      ce=.7
      eps=.1
cccccccccccc
      ngtp=1
      rm0=0.5
      itrn=8*60*60
      rm1=0.25
      itrn1=4*24*60*60
c      rm2=0.5
      rm2=0.25
      itrn2=4*24*60*60+12*60*60
      rm3=0.25
      cotor=1.e-5
      if (icps_bm .eq. 1) then
         rm0=0.5
         rm1=0.5
         rm2=0.5
         rm3=0.5
      endif

ccccc  total simulation time (min)

      iend=time*60.+.1
c
      d2t=dt
      sec=dt
c
      aminut=sec/60.
c     *****************
      ir36=21600
      r36=1./float(ir36)
c     *****************
      id=0
      idq=0
      n=1
      m=0
      m1=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc   base is the routine for initial conditions     ccccccccccccccccc
      call base (irs)
ctao
c      NUM_DIF=0
c      IF (NUM_DIF .EQ. 1) THEN
c        A2K=0.4E6
c        AAK=0.00175/DT
c      ENDIF
c
c     kdist=4
      kdist=2  ! cccshie by tao 5/3/01
c      if(.01*z2(3).ge.300.) kdist=3 ! cccshie by tao 5/3/01
c     print*,'kdist is',kdist        ! cccshie by tao 5/3/01
        kdistm=kdist
      if (icerh .eq. 1) call consat (rho)
      if (icelin .eq. 1) call consatl (rho)
c
      call slvpi (0)
c
      IF (ilif .EQ. 1) THEN
         rlat=0.0
         hrl=4.
         iday=18
         month=12
         kdist=3
         ice913=1
         iopcloud=1
         idom=1
         incar=0
      ENDIF

      iflag=0

      if (irad .eq. 1) then
         call zang (month,iday,cosz,rlat,hrl)
         call radrat (iflag,cosz,npp1,pl,pa)

C ------------- GCSS ---------------------------
    
        tmp1=0
        tmp2=0
        tmp3=0
        tmp4=0
        tmp5=0
        tmp6=0
        tmp7=0
        icount=0
c        print*,'printing out rflux'
            do i=2,nx/2-1
              icount=icount+1
              tmp1=rflux(i,1)+tmp1
              tmp2=rflux(i,2)+tmp2
              tmp3=rflux(i,3)+tmp3
              tmp4=rflux(i,4)+tmp4
              tmp5=rflux(i,5)+tmp5
              tmp6=rflux(i,6)+tmp6
              tmp7=rflux(i,7)+tmp7
            enddo
            tmp1=tmp1/real(icount)
            tmp2=tmp2/real(icount)
            tmp3=tmp3/real(icount)
            tmp4=tmp4/real(icount)
            tmp5=tmp5/real(icount)
            tmp6=tmp6/real(icount)
            tmp7=tmp7/real(icount)

c          print*,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7

C ----------------------------------------------------------------------

        iflag=1

      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (ISFC .EQ. 1 .AND. ITOGA .NE. 2) CALL PBLIN (0,isfcave)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if (iconflux .eq. 1) then
c           a1=1.e-3*rho1(2)
c          vvss=max(abs(ub(2)), 400.)
c          y111=a1*(1.1+.04*(.01*vvss))*vvss
c         y222=max(y111*(ts(2,1)-tb(2)),0.)
c         y333=max(y111*(qs(2,1)-qb(2)),0.)
c      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ndt_stat=60
       itime_ave=60*60
       mmave=1
c      if (irs .eq. 0) then
c        call dom_ave (aminut,tave,denq,thee)
c          avett(mmave)=tave
c          avetq(mmave)=denq
c          avet950(mmave)=thee
c          t_ub(mmave)=ub(2)
c          t_sfcq(mmave)=y333
c          t_sfct(mmave)=y222
c          t_microq(mmave)=0.
c          t_largeq(mmave)=0.
c          t_larget(mmave)=0.
c          t_microt(mmave)=0.
c      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       kt1=2
       kt2=kles
       it1=2
       it2=iles
       kc1=2
       kc2=kles
       ki1=2
       ki2=kles
       kr1=2
       kr2=kles
       ks1=2
       ks2=kles
       kg1=2
       kg2=kles

ccccc   adjustment for terminal velocity of rain
      do 106 k=2,kles
        fv(k)=sqrt(rho(2)*rrho(k))
        fv1(k)=sqrt(rho1(2)*rrho1(k))
  106 continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctao 9/22/97
c      this nuding coefficient should be used for NCAR model set-ups
c      imlifting = 0 and iuvbar = 0
c
       do k=1,nz
         ubi(k)=0.0
         vbi(k)=0.0
         do ii=1,itt3
           ub_2h(k,ii)=0.
           vb_2h(k,ii)=0.
         enddo
       enddo
       do ii=1,itt
          ip=ii+1
           if(ii .eq. itt) ip=itt
          iip3=(ii-1)*3+1
CHECK
          iipp3=iip3+1
          iippp3=iipp3+1
CHECK
         do k=1,nz
             y1(k)=(ub_6h(k,ip)-ub_6h(k,ii))/3.
             y2(k)=(vb_6h(k,ip)-vb_6h(k,ii))/3.
           ub_2h(k,iip3)=ub_6h(k,ii)
           ub_2h(k,iipp3)=ub_6h(k,ii)+y1(k)
           ub_2h(k,iippp3)=ub_6h(k,ii)+2.*y1(k)
           vb_2h(k,iip3)=vb_6h(k,ii)
           vb_2h(k,iipp3)=vb_6h(k,ii)+y2(k)
           vb_2h(k,iippp3)=vb_6h(k,ii)+2.*y2(k)
         enddo
       enddo
ctao (12-23-98)
        do k=1,nz
          ubi(k)=ub(k)
          vbi(k)=vb(k)
c         ubi(k)=ub_2h(k,2)
c         vbi(k)=vb_2h(k,2)
        enddo
ctao (12-23-98)
       do ii=1,itt
          iip3=(ii-1)*3+1
          iipp3=iip3+1
          iippp3=iipp3+1
         print*
         print*,'itt =',ii,iip3,iipp3,iippp3
         print*
         write(6,24100)
        do k1=2,kles
          k=kles+2-k1
          km=k-1
          write(6,20300) km,ub_2h(k,iip3),ub_2h(k,iipp3),ub_2h(k,iippp3)
     1                  ,vb_2h(k,iip3),vb_2h(k,iipp3),vb_2h(k,iippp3)
         enddo
       enddo

20300 format(i4,10f12.5,/)
24100 format(//,'level          u           u1          u2          v 
     1       v1          v2')

c      istrt_2h=0
      ir2hour=3600
      nudge=0
      r2hour=1./float(ir2hour)
      factor_nuding=r2hour
c
ctao 9/22/97
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       do k=1,nz
         tls(k)=0.0
         qls(k)=0.0
         UBT(K)=0.
         VBT(K)=0.
         WBT(K)=0.
         Q1T(K)=0.
         Q2T(K)=0.
         Q1ZT(K)=0.
         Q2ZT(K)=0.
         aq1t(k)=q1_6h(k,1)
         aq2t(k)=q2_6h(k,1)
         aq1zt(k)=q1z_6h(k,1)
         aq2zt(k)=q2z_6h(k,1)
       enddo
       IF (ilif .EQ. 0) Then
         do k=1,nz
           UBT(K)=(UB_6H(K,2)-UB_6H(K,1))*R36
           VBT(K)=(VB_6H(K,2)-VB_6H(K,1))*R36
           WBT(K)=(WB_6H(K,2)-WB_6H(K,1))*R36
           Q1T(K)=(Q1_6H(K,2)-Q1_6H(K,1))*R36
           Q2T(K)=(Q2_6H(K,2)-Q2_6H(K,1))*R36
           Q1ZT(K)=(Q1Z_6H(K,2)-Q1Z_6H(K,1))*R36
           Q2ZT(K)=(Q2Z_6H(K,2)-Q2Z_6H(K,1))*R36
c           aq1t(k)=0.0
c           aq2t(k)=0.0
c           aq1zt(k)=0.0
c           aq2zt(k)=0.0
           if (incar .eq. 1) then
             Q1ZT(K)=0.0
             Q2ZT(K)=0.0
           endif
         enddo
       ENDIF
         tsdt=(ts(1,2)-ts(1,1))*R36
         qsdt=(qs(1,2)-qs(1,1))*R36
         thsdt=(ths(1,2)-ths(1,1))*R36
         psdt=(pss(1,2)-pss(1,1))*R36
         do i=1,imax
           tairsfc(i)=ts(1,1)
           qairsfc(i)=qs(1,1)
           thairsf(i)=ths(1,1)
           pairsfc(i)=pss(1,1)
         enddo
      write(6,155)
      write(6,154) nrun,il2,kl2,dx,dz,time,dt,psfc,rd1,rd2,bound,
     1  ck,ce,rm0,rm1,eps,imp,iwp,imx,ised1,itrn,ngtp,cotor
c     ******************************************************************
      nran=0
      rm=rm0/150.
      do 107 i=2,iles
 107    rx(i)= (ran3(ised1)-0.5)*rm
      nran=nran+1
c     ***********************
      if (irs.eq.0 .and. n.eq.1) go to 6666
      if (irs.eq.0 .and. n.eq.0) go to 5555
       call wtap (dpt,7)
       call wtap (dqv,7)
       call wtap (qcl,7)
       call wtap (qrn,7)
       call wtap (qci,7)
       call wtap (qcs,7)
       call wtap (qcg,7)
       call wtap (u,7)
       call wtap (v,7)
       call wtap (w,7)
       call wtap (ak,7)
       call wtap (dpt1,7)
       call wtap (dqv1,7)
       call wtap (qcl1,7)
       call wtap (qrn1,7)
       call wtap (qci1,7)
       call wtap (qcs1,7)
       call wtap (qcg1,7)
       call wtap (uu1,7)
       call wtap (vv1,7)
       call wtap (ww1,7)
       call wtap (ak1,7)
       call wtap (umd,7)
       call wtap (vmd,7)
       call wtap (wmd,7)
       call wtap (rsw,7)
       call wtap (rlw,7)
       call wtap (ppress,7)
       call wtapt (7)
       call wstart (7)
c      call wstrad (7,cosz,month,iday,hrl)
       call wstrad (7,cosz,month,iday,hrl,rlat)  ! 8/21/01 shie
       call wstrad2 (7)  ! newly added! 8/21/01 shie
      dt=rdt
      d2t=dt+dt
      iend=time*60.+.1
      write(6,12345) n,isec,nran,sec,aminut,rdt

      do k=1,nz
c       wb(k)=0.0
c        taa(k)=qa(k)
      end do

      rm=rm0/150.
       if(isec.ge.itrn) rm=rm1/150.
       if(isec.ge.itrn1) rm=rm2/150.
       if(isec.ge.itrn2) rm=rm3/150.
      do 111 irk=1,nran-1
      do 111 i=2,iles
 111  rx(i)= (ran3(ised1)-0.5)*rm
      CALL SEPCA (ISEC,ICS,IV,ACO,ACO1,AAN,AAN1,LCONV,LANVL,LNSPT)
c     *****************************************************************
 3333 n=n+1

c      rbud = rbud + 1.

      call timing(time2)
      dtime=time2-time1
      time1=time2
      if(isec.eq.isec/60*60)then
        PRINT 241, N-1,SEC,AMINUT,DTIME
  241   FORMAT('TIMESTEP',I5,' MODEL TIME(S,MIN)=',F10.2,F10.2,3X,'CPU 
     1TIME(S)=',F8.3)
      endif
      sec=sec+dt
      aminut=sec/60.
      isec=sec+.1
      id=0
        if(isec.eq.isec/ndt_stat*ndt_stat) id=1
      idpbl=0
        if(isec.eq.isec/IPBLDT*IPBLDT) idpbl=1
        if (isfc .eq. 0) idpbl=0
        if (itoga .NE. 1) idpbl=0
      hrl=hrl+dt/3600.
      iradd=0
      if (irad .eq. 1) then
        if(isec.eq.isec/IPBLDT*IPBLDT) then
          call zang (month,iday,cosz,rlat,hrl)
          iradd=1
        endif
      endif
c
cc
ctao
c
cc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (ilif .EQ. 1) THEN
c       fix sst=29.0
	tsfcLIF=29.0
	tslLIF=(tsfcLIF+273.16)/pi(2)
	qslLIF=3.799052e3/p0(2)*exp(17.26939*tsfcLIF/(tsfcLIF+237.3))
        qslLIF=qslLIF*.975
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (icase .eq. 0) then
        ICY_UB=ISEC/ir36+1
        DO ICY=1,ICY_UB
          ICHK=ISTRT+(ICY-1)*ir36
          IF(ISEC.EQ.ICHK)THEN
            ITFOR=ICY
            IF (ilif .EQ. 0)ITFOR=MIN(ICY+1,ITT-1)
            IF (IWBAR .EQ. 0) THEN
              DO K=2,KLES
                UBT(K)=(UB_6H(K,itfor+1)-UB_6H(K,itfor))*R36
                VBT(K)=(VB_6H(K,itfor+1)-VB_6H(K,itfor))*R36
                WBT(K)=(WB_6H(K,itfor+1)-WB_6H(K,itfor))*R36
                Q1T(K)=(Q1_6H(K,itfor+1)-Q1_6H(K,itfor))*R36
                Q2T(K)=(Q2_6H(K,itfor+1)-Q2_6H(K,itfor))*R36
                Q1ZT(K)=(Q1Z_6H(K,itfor+1)-Q1Z_6H(K,itfor))*R36
                Q2ZT(K)=(Q2Z_6H(K,itfor+1)-Q2Z_6H(K,itfor))*R36
                 aq1t(k)=q1_6h(k,itfor)
                 aq2t(k)=q2_6h(k,itfor)
                 aq1zt(k)=q1z_6h(k,itfor)
                 aq2zt(k)=q2z_6h(k,itfor)
              ENDDO
              tsdt=(ts(1,itfor+1)-ts(1,itfor))*R36
              qsdt=(qs(1,itfor+1)-qs(1,itfor))*R36
              thsdt=(ths(1,itfor+1)-ths(1,itfor))*R36
              psdt=(pss(1,itfor+1)-pss(1,itfor))*R36
               do i=1,imax
                 tairsfc(i)=ts(1,itfor)
                 qairsfc(i)=qs(1,itfor)
                 thairsf(i)=ths(1,itfor)
                 pairsfc(i)=pss(1,itfor)
               enddo
            ELSE
              DO K=2,KLES
                if (iuvbar .eq. 0) UBT(K)=0.
                if (iuvbar .eq. 0) VBT(K)=0.
                WBT(K)=0.
                Q1T(K)=0.
                Q2T(K)=0.
              ENDDO
              tsdt=0.
              qsdt=0.
              thsdt=0.
              psdt=0.
            ENDIF
            WRITE (6,*) 'UPDATE UBT AT ', ISEC,' ICY=',itfor
            WRITE (6,*) 'UPDATE WBT AT ', ISEC,' ICY=',itfor
          ENDIF
        ENDDO
      endif
C-----------------------------------------------------------------------
      do 631 i=2,iles
       y1(i)=(dpt1(i,1)+ta1(1))*(1.+.61*(dqv1(i,1)+qa1(1))
     1        -qcl1(i,1)-qrn1(i,1)-qci1(i,1)-qcs1(i,1)-qcg1(i,1))
  631  y2(i)=(dpt1(i,2)+ta1(2))*(1.+.61*(dqv1(i,2)+qa1(2))
     1        -qcl1(i,2)-qrn1(i,2)-qci1(i,2)-qcs1(i,2)-qcg1(i,2))
      do 632 k=2,kles
        kp=k+1
      do 632 i=2,iles
       y3(i)=(dpt1(i,kp)+ta1(kp))*(1.+.61*(dqv1(i,kp)+qa1(kp))
     1    -qcl1(i,kp)-qrn1(i,kp)-qci1(i,kp)-qcs1(i,kp)-qcg1(i,kp))
       if(y3(i)-y1(i) .ge. 0.) then
         xxw(i,k)=0.0
       else
         xxw(i,k)=0.5
       endif
       y1(i)=y2(i)
       y2(i)=y3(i)
  632 continue
      do 634 i=2,iles
       xxw(i,1)=xxw(i,2)
  634  xxw(i,kmax)=xxw(i,kles)
      do 635 k=2,kmax
       y1(k)=f0(k)+f0(k-1)
  635  y2(k)=.5*y1(k)
      do 636 k=3,kmax
      do 636 i=2,iles
       df0(i,k)=2.17375*y1(k)*(dqv(i,k)+dqv(i,k-1)+qa(k)+qa(k-1))
     1          /(dpt(i,k)+dpt(i,k-1)+ta(k)+ta(k-1))**2
       dfc(i,k)=-am1(k)*(ak(i,k)*(1.+xxw(i,k))+ak(i,k-1)
     1                          *(1.+xxw(i,k-1)))
     2        *(dpt1(i,k)-dpt1(i,k-1)+ta1(k)-ta1(k-1)+y2(k)
     3        *(dqv1(i,k)-dqv1(i,k-1)+qa1(k)-qa1(k-1)))
     4        /((1.+df0(i,k)*y2(k))*dz)
  636 continue
      do 637 i=2,iles
       df0(i,2)=0.
  637  dfc(i,2)=0.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (IDPBL .EQ. 1) THEN
           do i=2,iles
              if (iricont .eq. 0) then
                ri180(i)=0.
              else
                ri180(i)=ri180(i)/float(iricont)
              endif
              riold180(i) = ri180(i)
           enddo
         CALL PBLIN (1,isfcave)
           iricont=0
           do i=2,iles
             ri180(i)=0.
           enddo
      ENDIF
c     *******************************************************
       it1=2
       it2=iles
       kt1=2
       kt2=kles
       irfg=0
       idq=0
       irf=1
       ismg=1
       IADVH=IADVH1
        IF (IJKADV .EQ. 1)  IADVH=2
      call advect (dpt,dpt1,ta1,TA,id,id, 0, 0, 0,id, 1., 1.)
      do 639 k=3,kmax
      do 639 i=2,iles
  639  dfc(i,k)=dfc(i,k)*df0(i,k)
      irfg=0
      irf=1
c     *******************************************************
      ismg=2
      call advect (dqv,dqv1,qa1,QA, 0, 0,id, 0, 0, 0, 1., 1.)
      irf=0
      do 640 k=3,kmax
      do 640 i=2,iles
       dfc(i,k)=-dfc(i,k)-am1(k)*(ak(i,k)+ak(i,k-1))
     1  *(dqv1(i,k)-dqv1(i,k-1)+qa1(k)-qa1(k-1)+qcl1(i,k)
     2  -qcl1(i,k-1)+qci1(i,k)-qci1(i,k-1))*rdz
  640 continue
c     *******************************************************
        it1=2
        it2=iles
        kc1=2
        kc2=kles
        ki1=2
        ki2=kles
        kr1=2
        kr2=kles
        ks1=2
        ks2=kles
        kg1=2
        kg2=kles
       if (idom .eq. 1) then
        do 305 k=2,kles
        do 305 i=2,iles
  305     xxx(i,k)=qcl(i,k)+qcl1(i,k)
        call domain (xxx,it1,it2,kc1,kc2)
        do 325 k=2,kles
        do 325 i=2,iles
          xxx(i,k)=qrn(i,k)+qrn1(i,k)
 325    continue
        call domain (xxx,it1,it2,kr1,kr2)
        IF (ICE .EQ. 1) THEN
          do 315 k=2,kles
          do 315 i=2,iles
            xxx(i,k)=qci(i,k)+qci1(i,k)
 315      continue
          call domain (xxx,it1,it2,ki1,ki2)
          do 335 k=2,kles
          do 335 i=2,iles
            xxx(i,k)=qcs(i,k)+qcs1(i,k)
 335      continue
          call domain (xxx,it1,it2,ks1,ks2)
          do 345 k=2,kles
          do 345 i=2,iles
            xxx(i,k)=qcg(i,k)+qcg1(i,k)
 345      continue
        call domain (xxx,it1,it2,kg1,kg2)
        ENDIF
       endif
c     **********************************************************
      idq=id
       kt1=kc1
       kt2=kc2
       ismg=3
       irf=0
       do 355 k=1,kmax
       do 355 i=2,iles
  355   qcl_buffer(i,k)=qcl(i,k)
      call advect (qcl,qcl1,qaA1,QAA1, 0, 0, 0, 1, 0, 0, 1., 0.)
c     **********************************************************
      IF (ICE .EQ. 1) THEN
        kt1=ki1
        kt2=ki2
        irf=0
        ismg=4
        call advect (qci,qci1,qAa1,QAA1, 0, 0, 0, 1, 0, 0, 0., 1.)
      ENDIF
c     **********************************************************
      do 641 k=1,kmax
      do 641 i=1,imax
  641  xxx(i,k)=ww1(i,k)
c
       irsg=0
      if (icerh .eq. 1) call terv (irsg,rho1,fv1)
      if (icelin .eq. 1) call tervl (irsg,rho1,fv1)
c
      do 701 i=2,iles
cc     RI(I)=36.E3*RHO1(2)*QRN(I,2)*WW1(I,2)
       RI(I)=36.E3*RHO1(2)*QRN(I,2)*VTP(I,2)
       ri180(i)=ri180(i)+ri(i)
       ri120ar(i)=ri120ar(i)+ri(i)
        if (icps_bm .eq. 1)  ri180(i)=ri180(i)+rainco(i,1)
       ar_c(i)=ar_c(i)+rainco(i,1)*dt*2.77778e-4
  701  ar(i)=ar(i)+ri(i)*dt*2.77778e-4
       IRICONT=IRICONT+1
       IRICONT1=IRICONT1+1
c
         if (ijkadv .eq. 1) then
           do 702 k=1,kmax
           do 702 i=1,imax
             ww1(i,k)=.5*(3.*xxx(i,k)-wmd(i,k))
  702      continue
         else
           do k=1,kmax
           do i=1,imax
             ww1(i,k)=xxx(i,k)-vtp(i,k)
           enddo
           enddo
         endif

       kt1=kr1
       kt2=kr2
      ismg=5
      call advect (qrn,qrn1,qaA1,QAA1, 0, 0, 0, 1, 1, 0, 0., 0.)
c     **********************************************************
      IF (ICE .EQ. 1) THEN
        irsg=1
        if (icerh .eq. 1) call terv (irsg,rho1,fv1)
        if (icelin .eq. 1) call tervl (irsg,rho1,fv1)
         if (ijkadv .eq. 1) then
           do 703 k=1,kmax
           do 703 i=1,imax
             ww1(i,k)=.5*(3.*xxx(i,k)-wmd(i,k))
  703      continue
         else
           do k=1,kmax
           do i=1,imax
             ww1(i,k)=xxx(i,k)-vtp(i,k)
           enddo
           enddo
         endif
        kt1=ks1
        kt2=ks2
        call advect (qcs,qcs1,qaA1,QAA1, 0, 0, 0, 1, 1, 0, 0., 0.)
c     **********************************************************
        irsg=2
        if (icerh .eq. 1) call terv (irsg,rho1,fv1)
        if (icelin .eq. 1) call tervl (irsg,rho1,fv1)
         if (ijkadv .eq. 1) then
           do 704 k=1,kmax
           do 704 i=1,imax
             ww1(i,k)=.5*(3.*xxx(i,k)-wmd(i,k))
  704      continue
         else
           do k=1,kmax
           do i=1,imax
             ww1(i,k)=xxx(i,k)-vtp(i,k)
           enddo
           enddo
         endif
        kt1=kg1
        kt2=kg2
        call advect (qcg,qcg1,qaA1,QAA1, 0, 0, 0, 1, 1, 0, 0., 0.)

      END IF
c     **********************************************************
      irf=1
      smr=0.
      smrc=0.
      do 73 i=2,iles
         smr=smr+ar(i)
         smrc=smrc+ar_c(i)
  73  continue
      smr=smr*ril2
      smrc=smrc*ril2

      do 705 k=1,kmax
      do 705 i=1,imax
  705  ww1(i,k)=xxx(i,k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       D22T=D2T
       IF(IJKADV .EQ. 1) D22T=DT
        afact=1.
        if (ijkadv .eq. 0) afact=0.5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       do k=2,kles
         y1(k)=d22t*rdz*rrho(k)
         y2(k)=.5*d22t*fcorg*rdz*am(k)*tb(k)
         IF (ICASE .EQ. 1) THEN
           Q1ZT(K)=tlsw(k)*1.004*3600.0*pi(k)
           Q2ZT(K)=qlsw(k)*2500.0*3600.0
           Q1T(K)=0.
           Q2T(K)=0.
           if (isui .eq. 1) Q1ZT(K)=0.
           if (isui .eq. 1) Q2ZT(K)=0.
         ENDIF
       enddo
       rtfor=d22t/(24.*3600.)
       rqfor=d22t/(24.*3600.)
        if (incar .eq. 0) then
          rqfor=d22t/(2500.*24.*3600.)
        endif
      do 708 k=2,kles
        kp=k+1
        km=k-1
         a3=0.
         a4=0.
         a5=0.
         a6=0.
         a7=0.
         a8=0.

        if(ijkadv .eq. 1) then
          do i=2,iles
             y3(i)=.5*(3.*ww1(i,k)-wmd(i,k))
             y4(i)=.5*(3.*ww1(i,kp)-wmd(i,kp))
          enddo
        else
          do i=2,iles
             y3(i)=ww1(i,k)
             y4(i)=ww1(i,kp)
          enddo
        endif
        DO 7071 I=2,ILES

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (ISNGTAO .EQ. 0) THEN

c           q1_g_v(i,k)=q1_g_v(i,k)-(AM1(K)*FD(K)*Y3(I)
c     1                  +AM1(KP)*FD(KP)*Y4(I))*Y1(K)
c
c           q2_g_v(i,k)=q2_g_v(i,k)-(AM1(K)*FE(K)*Y3(I)
c     1                  +AM1(KP)*FE(KP)*Y4(I))*Y1(K)

         DPT(I,K)=DPT(I,K)-(AM1(K)*FD(K)*Y3(I)
     1            +AM1(KP)*FD(KP)*Y4(I))*Y1(K)
c    2            +Y2(K)*(VV1(I,K)-VB1(K))*(UU1(I,KP)-UU1(I,KM))
     3            -y0(i)*aQ1T(K)*rtfor/PI(K)
c    4            -y0(i)*aQ1ZT(K)*rtfor/PI(K)
         DQV(I,K)=DQV(I,K)-(AM1(K)*FE(K)*Y3(I)
     1            +AM1(KP)*FE(KP)*Y4(I))*Y1(K)
     2            -y0(i)*aQ2T(K)*rqfor
c    3            -y0(i)*aQ2ZT(K)*rqfor

          a1=-(am1(k)*fd(k)*y3(i)+am1(kp)*fd(kp)*y4(i))*y1(k)
     1       -y0(i)*aq1t(k)*rtfor/pi(k)-y0(i)*aq1zt(k)*rtfor/pi(k)
          a2=-(am1(k)*fe(k)*y3(i)+am1(kp)*fe(kp)*y4(i))*y1(k)
     1       -y0(i)*aq2t(k)*rqfor-y0(i)*aq2zt(k)*rqfor

         a5=a5-(AM1(K)*FD(K)*Y3(I)+AM1(KP)*FD(KP)*Y4(I))*Y1(K)
     1        -y0(i)*aQ1T(K)*rtfor/PI(K)
         a6=a6-(AM1(K)*FE(K)*Y3(I)+AM1(KP)*FE(KP)*Y4(I))*Y1(K)
     2        -y0(i)*aQ2T(K)*rqfor

          a77=-y0(i)*aq1zt(k)*rtfor/pi(k)
          a88=-y0(i)*aq2zt(k)*rqfor
       ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dpt(i,k)=dpt(i,k)-(am1(k)*fd(k)*(Y3(I)+wb(k)*y0(i))
     1            +am1(kp)*fd(kp)*(Y4(I)+wb(kp)*y0(i)))*y1(k)
c    2            +y2(k)*(vv1(i,k)-vb1(k))*(uu1(i,kp)-uu1(i,km))
     3            -y0(i)*aQ1T(K)*rtfor/PI(K)
         dqv(i,k)=dqv(i,k)-(am1(k)*fe(k)*(Y3(I)+wb(k)*y0(i))
     1            +am1(kp)*fe(kp)*(Y4(I)+wb(kp)*y0(i)))*y1(k)
     2            -y0(i)*aQ2T(K)*rqfor
          a1=-(am1(k)*fd(k)*(Y3(I)+wb(k)*y0(i))
     1        +am1(kp)*fd(kp)*(Y4(I)+wb(kp)*y0(i)))*y1(k)
     2       -y0(i)*aQ1T(K)*rtfor/PI(K)
          a2=-(am1(k)*fe(k)*(Y3(I)+wb(k)*y0(i))
     1        +am1(kp)*fe(kp)*(Y4(I)+wb(kp)*y0(i)))*y1(k)
     2       -y0(i)*aQ2T(K)*rqfor
c
          a5=a5-(am1(k)*fd(k)*(Y3(I)+wb(k)*y0(i))
     1          +am1(kp)*fd(kp)*(Y4(I)+wb(kp)*y0(i)))*y1(k)
     2         -y0(i)*aQ1T(K)*rtfor/PI(K)
          a6=a6-(am1(k)*fe(k)*(Y3(I)+wb(k)*y0(i))
     1        +am1(kp)*fe(kp)*(Y4(I)+wb(kp)*y0(i)))*y1(k)
     2       -y0(i)*aQ2T(K)*rqfor
c
          a77=-(am1(k)*fd(k)*wb(k)*y0(i)
     1         +am1(kp)*fd(kp)*wb(kp)*y0(i))*y1(k)
          a88=-(am1(k)*fe(k)*wb(k)*y0(i)
     1         +am1(kp)*fe(kp)*wb(kp)*y0(i))*y1(k)

        ENDIF
c         a5=a5-(am1(k)*fd(k)*y3(i)+am1(kp)*fd(kp)*y4(i))*y1(k)
c         a6=a6-(am1(k)*fe(k)*y3(i)+am1(kp)*fe(kp)*y4(i))*y1(k)
         a7=a7+a77
         a8=a8+a88
         a3=a3+a1
         a4=a4+a2
 7071   CONTINUE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ttlsw(k)=ttlsw(k)+a5*ril2*afact
       tqlsw(k)=tqlsw(k)+a6*ril2*afact
       ttlsw1(k)=ttlsw1(k)+a7*ril2*afact
       tqlsw1(k)=tqlsw1(k)+a8*ril2*afact
c       tls(k)=tls(k)+.5*a3*ril2
c       qls(k)=qls(k)+.5*a4*ril2
       tls(k)=tls(k)+a3*ril2*afact
       qls(k)=qls(k)+a4*ril2*afact
  708 continue
c     ********************************
      if(isec.ne.isec/300*300) go to 6666
      if(isec.ge.itrn) rm=rm1/150.
      if(isec.ge.itrn1) rm=rm2/150.
      if(isec.ge.itrn2) rm=rm3/150.
      do 16 i=2,iles
        y1(i)=rx(i)
        rx(i)= (ran3(ised1)-0.5)*rm
   16  continue
c
       IF (IJKADV .EQ. 0) THEN
          do k=kdistm,kdist
          do i=2,iles
             DPT(I,k)=DPT(I,k)+(Y1(I)-RX(I))*DT
          enddo
          enddo
       ENDIF
c
      nran=nran+1
c     ***************************************
 6666 continue
        D22T=D2T
         IF(IJKADV .EQ. 1) D22T=DT
      do k=kdistm,kdist
      DO 18 I=2,ILES
   18   DPT(I,k)=DPT(I,k)+RX(I)*D22T
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       if (isec .ge. 3*3600) then
       if (isfc .eq. 1 .AND. ITOGA .EQ. 2) then
        D22T=D2T
         IF(IJKADV .EQ. 1) D22T=DT
       A1=1.0E-3*RHO1(2)
       A2=D22T*AM(2)*RDZ*RRHO(2)
       A3=0.
       A4=0.
       DO 300 I=2,ILES
        IP=I+1
         IF (I .EQ. ILES) IP=2
        IF (isfcave .EQ. 0) THEN
          vvss=sqrt(.25*(uu1(i,2)+uu1(ip,2))**2 
     1            + .25*(vv1(i,2)+vv1(ip,2))**2)
          if(ijkadv .eq. 1) vvss=sqrt((.25*(3.*uu1(i,2)-umd(i,2)
     1         +3.*uu1(ip,2)-umd(ip,2)))**2
     2                              + (.25*(3.*vv1(i,2)-vmd(i,2)
     3         +3.*vv1(ip,2)-vmd(ip,2)))**2)
        else
          vvss=sqrt( ub1(2)**2 + vb1(2)**2 )
        endif
        IF (ilif .EQ. 1) VVSS=abs( u(i,2) )

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         Y4(I)=MAX(VVSS, 100.0)
          Y1(I)=A1*(1.1+.04*(.01*Y4(I)))*Y4(I)
        IF (isfcave .EQ. 0) THEN
          Y2(I)=MAX(Y1(I)*(thairsf(I)-TA1(2)-DPT(I,2)) ,0.)
          Y3(I)=MAX(Y1(I)*(fact_qvs*qairsfc(I)-QA1(2)-DQV(I,2)) ,0.)
        ELSE
          Y2(I)=MAX(Y1(I)*(thairsf(I)-TA1(2)) ,0.)
          Y3(I)=MAX(Y1(I)*(fact_qvs*qairsfc(I)-QA1(2)) ,0.)
        ENDIF
        IF (ilif .EQ. 1) THEN
          Y2(I)=MAX(Y1(I)*(tslLIF-TA1(2)-DPT(I,2)) ,0.)
          Y3(I)=MAX(Y1(I)*(qslLIF-QA1(2)-DQV(I,2)) ,0.)
        ENDIF
C SAVE SENSIBLE AND LATENT HEAT FLUXES
        SUN_4(I,3)=Y2(I)
        SUN_4(I,4)=Y3(I)
         tsfc_1(i)=y2(i)
         qsfc_1(i)=y3(i)
c***********************************************************************

c         q1_d_v(i,2)=q1_d_v(i,2)+a2*y2(i)
c         q2_d_v(i,2)=q2_d_v(i,2)+a2*y3(i)

        dpt(i,2)=dpt(i,2)+a2*y2(i)
        dqv(i,2)=dqv(i,2)+a2*y3(i)
        a3=a3+y2(i)
  300   a4=a4+y3(i)
      stf(2)=stf(2)+a3*dt*ril2
      sqf(2)=sqf(2)+a4*dt*ril2
c      stfff2=A3*DT*RIL2
c      sqfff2=10.*A4*DT*RIL2
      endif
C      endif
cc   ***  churchill and houzes' convective and anvil separation   *****
      IF (ID .EQ. 1) THEN
        CALL SEPCA (ISEC,ICS,IV,ACO,ACO1,AAN,AAN1,LCONV,LANVL,LNSPT)
      ENDIF
cc    ***   insert radiative transfer processes here *****
      if(irad.eq.1) then
        iflag=1
        ircl=ircl+1
          if(isec .eq. isec/21600*21600) iflag=2
        if (iradd .eq. 1) then
          call radrat (iflag,cosz,npp1,pl,pa)
C ------------- GCSS ---------------------------
    
        tmp1=0
        tmp2=0
        tmp3=0
        tmp4=0
        tmp5=0
        tmp6=0
        tmp7=0
        icount=0

c        print*,'printing out rflux'
            do i=1,nx/2-1
              icount=icount+1
              tmp1=rflux(i,1)+tmp1
              tmp2=rflux(i,2)+tmp2
              tmp3=rflux(i,3)+tmp3
              tmp4=rflux(i,4)+tmp4
              tmp5=rflux(i,5)+tmp5
              tmp6=rflux(i,6)+tmp6
              tmp7=rflux(i,7)+tmp7
            enddo
            tmp1=tmp1/real(icount)
            tmp2=tmp2/real(icount)
            tmp3=tmp3/real(icount)
            tmp4=tmp4/real(icount)
            tmp5=tmp5/real(icount)
            tmp6=tmp6/real(icount)
            tmp7=tmp7/real(icount)

c          print*,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
         endif

C ---------------------------------------------------------------
        D22T=D2T
         IF(IJKADV .EQ. 1) D22T=DT
        do 31 k=2,kles
        do 31 i=2,iles
c          q1_rad(i,k)=q1_rad(i,k)+(RSW(I,K)-RLW(I,K))*d22t
   31     dpt(i,k)=dpt(i,k)+(rsw(i,k)-rlw(i,k))*d22t
      endif
cc    **************************
        it1=2
        it2=iles
        kt1=2
        kt2=kles
      if (idom .eq. 1) then
        t0=273.16
        do 330 k=2,kles
          rp0=3.799052e3/p0(k)
        do 330 i=2,iles
         tair(i)=(dpt(i,k)+ta1(k))*pi(k)
         if(tair(i) .ge. t0) then
           qss(i)=rp0*exp(c172-c409/(t0-c358))
         else
           qss(i)=rp0*exp(c218-c580/(tair(i)-c76))
         endif
         IF (ICE .EQ. 0) qss(i)=rp0*exp(c172-c409/(t0-c358))
         qv(i)=max(0.E0, dqv(i,k)+qa1(k)-qss(i))
          if (qcl(i,k).le.1.e-20) qcl(i,k)=0.0
          if (qrn(i,k).le.1.e-20) qrn(i,k)=0.0
          if (qci(i,k).le.1.e-20) qci(i,k)=0.0
          if (qcs(i,k).le.1.e-20) qcs(i,k)=0.0
          if (qcg(i,k).le.1.e-20) qcg(i,k)=0.0
  330   xxx(i,k)=qv(i)+qcl(i,k)+qci(i,k)+qrn(i,k)+qcs(i,k)+qcg(i,k)
        call domain (xxx,it1,it2,kt1,kt2)
      endif
c    **************************
      IF (ICE .EQ. 0) call sat (n,id,fv)
      IF (ICE .EQ. 1 .and. icerh .eq. 1) call satice(n,id,KCPS,fv)
      IF (ICE .EQ. 1 .and. icelin .eq. 1) CALL SATICEL (n,id,fv)
C *******************************************
      do k=2,kles
        aq1t(k)=aq1t(k)+q1t(k)*dt
        aq2t(k)=aq2t(k)+q2t(k)*dt
        aq1zt(k)=aq1zt(k)+q1zt(k)*dt
        aq2zt(k)=aq2zt(k)+q2zt(k)*dt
        if (incar .eq. 1) then
          aq1zt(k)=0.0
          aq2zt(k)=0.0
        endif
        wb(k)=wb(k)+wbt(k)*dt
ctao (12-23-98)
        ubi(k)=ubi(k)+ubt(k)*dt
        vbi(k)=vbi(k)+vbt(k)*dt
ctao (12-23-98)
      enddo

      do i=1,imax
        tairsfc(i)=tairsfc(i)+tsdt*dt
        qairsfc(i)=qairsfc(i)+qsdt*dt
        thairsf(i)=thairsf(i)+thsdt*dt
        pairsfc(i)=pairsfc(i)+psdt*dt
      enddo
          iihalf=il2/2+1
         if(isec.eq. isec/10800*10800)
     1       print*,aminut,tairsfc(iihalf),qairsfc(iihalf)*1000.,
     2                pairsfc(iihalf)/1000.,thairsf(iihalf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctao 9/22/97
ctao 12/23/98
c       if (isec .eq. isec/ir2hour*ir2hour) then
c         istrt_2h=istrt_2h+1
c         do k=1,nz
c           ubi(k)=ub_2h(k,istrt_2h)
c           vbi(k)=vb_2h(k,istrt_2h)
c         enddo
c         print*,aminut,isec
c         do k=2,kles
c           print*,k,ubi(k),vbi(k)
c         enddo
c       endif
c
ctao 12/23/98
ctao 9/22/97
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       if(isec.eq.isec/ibudsec*ibudsec) then
c         do k=1,nz
c           do i=1,nx
c             q1a_g_h(i,k) = q1_g_h(i,k) / rbud
c             q1a_g_v(i,k) = q1_g_v(i,k) / rbud
c             q1a_d_h(i,k) = q1_d_h(i,k) / rbud
c             q1a_d_v(i,k) = q1_d_v(i,k) / rbud
c             q2a_g_h(i,k) = q2_g_h(i,k) / rbud
c             q2a_g_v(i,k) = q2_g_v(i,k) / rbud
c             q2a_d_h(i,k) = q2_d_h(i,k) / rbud
c             q2a_d_v(i,k) = q2_d_v(i,k) / rbud
c             q1a_rad(i,k) = q1_rad(i,k) / rbud
c
c             q1_g_h(i,k)=0.
c             q1_g_v(i,k)=0.
c             q1_d_h(i,k)=0.
c             q1_d_v(i,k)=0.
c             q2_g_h(i,k)=0.
c             q2_g_v(i,k)=0.
c             q2_d_h(i,k)=0.
c             q2_d_v(i,k)=0.
c             q1_rad(i,k)=0.
c
c           enddo
c         enddo
c         rbud = 0.
c       endif

cc    ******************************************************************
      d2t=dt+dt
      isec=sec+.1
      smt=0.
      smv=0.
      smq=0.
      smc=0.
      sme=0.
      smqc=0.
      smqr=0.
      smqi=0.
      smqs=0.
      smqg=0.
c
      sddd=0.
      ssss=0.
      shhh=0.
      sccc=0.
      smmm=0.
      sfff=0.
c
      do 32 k=2,kles
        y1(k)=rho(k)*dz*10./am(k)
        y2(k)=y1(k)/f0(k)
       smt=smt+st(k)*y2(k)
       smv=smv+sv(k)*y1(k)
       smq=smq+sq(k)*y1(k)
       smc=smc+sc(k)*y1(k)
       sme=sme+se(k)*y1(k)
c
       sddd=sddd+s_dep(k)*y1(k)
       ssss=ssss+s_sub(k)*y1(k)
       shhh=shhh+s_qrs(k)*y1(k)
       sccc=sccc+s_qrl(k)*y1(k)
       smmm=smmm+s_mel(k)*y1(k)
       sfff=sfff+s_frz(k)*y1(k)
c
       smqc=smqc+sqc(k)*y1(k)
       smqr=smqr+sqr(k)*y1(k)
       smqi=smqi+sqi(k)*y1(k)
       smqs=smqs+sqs(k)*y1(k)
   32  smqg=smqg+sqg(k)*y1(k)
c     *********************************
      if(isec.lt.iwtp) go to 1111
      if(isec.ne.isec/iwp*iwp) go to 1111

       call wtap2 (ak,bnxnz,nx,nz,4)
       call wtap2 (dpt,bnxnz,nx,nz,4)
       call wtap2 (dqv,bnxnz,nx,nz,4)
       call wtap2 (qcl,bnxnz,nx,nz,4)
       call wtap2 (qrn,bnxnz,nx,nz,4)
       call wtap2 (qci,bnxnz,nx,nz,4)
       call wtap2 (qcs,bnxnz,nx,nz,4)
       call wtap2 (qcg,bnxnz,nx,nz,4)
       call wtap2 (u,bnxnz,nx,nz,4)
       call wtap2 (v,bnxnz,nx,nz,4)
       call wtap2 (w,bnxnz,nx,nz,4)
       call wtap2 (rsw,bnxnz,nx,nz,4)
       call wtap2 (rlw,bnxnz,nx,nz,4)

       CALL WTAP_UV

C     TAO 12-10-97
        if(icps_bm .eq. 1) then
        CALL WTAP2 (TMODO,BNXNZ,NX,NZ,4)
        CALL WTAP2 (QMODO,BNXNZ,NX,NZ,4)
        endif
C     TAO 12-10-97

       call wtap1 (ub,bnz,nz,4)
       call wtap1 (vb,bnz,nz,4)
       call wtap1 (ub1,bnz,nz,4)
       call wtap1 (vb1,bnz,nz,4)
       call wtap1 (tb,bnz,nz,4)
       call wtap1 (qb,bnz,nz,4)
       call wtap1 (rho,bnz,nz,4)
       call wtap1 (rho1,bnz,nz,4)
       call wtap1 (ta,bnz,nz,4)
       call wtap1 (qa,bnz,nz,4)
       call wtap1 (ta1,bnz,nz,4)
       call wtap1 (qa1,bnz,nz,4)
       call wtap1 (pi,bnz,nz,4)
       call wtap1 (p0,bnz,nz,4)
       call wtap1 (wb,bnz,nz,4)
       call wtap1 (wbt,bnz,nz,4)
       call wtap1 (fd,bnz,nz,4)
       call wtap1 (fe,bnz,nz,4)
       call wtap1 (sqc,bnz,nz,4)
       call wtap1 (sqr,bnz,nz,4)
       call wtap1 (sqi,bnz,nz,4)
       call wtap1 (sqs,bnz,nz,4)
       call wtap1 (sqg,bnz,nz,4)
       call wtap1 (stqc,bnz,nz,4)
       call wtap1 (stqr,bnz,nz,4)
       call wtap1 (stqi,bnz,nz,4)
       call wtap1 (stqs,bnz,nz,4)
       call wtap1 (stqg,bnz,nz,4)
       call wtap1 (aq1t,bnz,nz,4)
       call wtap1 (aq2t,bnz,nz,4)
       call wtap1 (aq1zt,bnz,nz,4)
       call wtap1 (aq2zt,bnz,nz,4)
       call wtap1 (suw,bnx,nx,4)
       call wtap1 (svw,bnx,nx,4)
       call wtap1 (swt,bnx,nx,4)
       call wtap1 (swq,bnx,nx,4)
       call wtap1 (tairsfc,bnx,nx,4)
       call wtap1 (qairsfc,bnx,nx,4)
       call wtap1 (thairsf,bnx,nx,4)
       call wtap1 (pairsfc,bnx,nx,4)
       call wtap2 (sun_4,bnx4,nx,4,4)
       call wtap2 (rflux,bnx7,nxi/2,7,4)

c       call wtap2(q1a_g_h,bnxnz,nx,nz,4)
c       call wtap2(q1a_g_v,bnxnz,nx,nz,4)
c       call wtap2(q1a_d_h,bnxnz,nx,nz,4)
c       call wtap2(q1a_d_v,bnxnz,nx,nz,4)
c       call wtap2(q2a_g_h,bnxnz,nx,nz,4)
c       call wtap2(q2a_g_v,bnxnz,nx,nz,4)
c       call wtap2(q2a_d_h,bnxnz,nx,nz,4)
c       call wtap2(q2a_d_v,bnxnz,nx,nz,4)
c       call wtap2(q1a_hyd,bnxnz,nx,nz,4)
c       call wtap2(q2a_hyd,bnxnz,nx,nz,4)
c       call wtap2(q1a_rad,bnxnz,nx,nz,4)


       call wtap0 (smr,4)
       call wtap0 (smrc,4)
       call wtap1 (ri,bnx,nx,4)
       call wtap1 (ar_c,bnx,nx,4)
       call wtap1 (ar,bnx,nx,4)
       call wtap1 (ri180,bnx,nx,4)
       call wtap1 (ri120ar,bnx,nx,4)

    
cc    *********************************
 1111 if(isec.lt.iwtp) go to 2222
      if(isec.ne.isec/imp*imp) go to 2222
      call wtap3 (S1,b16nz14,16,nz,14,4)
      call wtap3 (SN1,b5nz14,5,nz,14,4)
      call wtap2 (S9,b16nz,16,nz,4)
      call wtap2 (SN9,b5nz,5,nz,4)
      call wtap2 (S10,b16nz,16,nz,4)
      call wtap2 (SN10,b5nz,5,nz,4)
      call wtap2 (S11,b16nz,16,nz,4)
      call wtap2 (SN11,b5nz,5,nz,4)
      call wtap2 (S12,b16nz,16,nz,4)
      call wtap2 (SN12,b5nz,5,nz,4)
      call wtap2 (S13,b16nz,16,nz,4)
      call wtap2 (SN13,b5nz,5,nz,4)
      call wtap2 (S14,b16nz,16,nz,4)
      call wtap2 (SN14,b5nz,5,nz,4)
      call wtap2 (S15,b16nz,16,nz,4)
      call wtap2 (SN15,b5nz,5,nz,4)
      call wtap2 (S16,b16nz,16,nz,4)
      call wtap2 (SN16,b5nz,5,nz,4)
      call wtap2 (S17,b16nz,16,nz,4)
      call wtap2 (SN17,b5nz,5,nz,4)
      call wtap2 (S18,b16nz,16,nz,4)
      call wtap2 (SN18,b5nz,5,nz,4)
      call wtap2 (S19,b16nz,16,nz,4)
      call wtap2 (SN19,b5nz,5,nz,4)
      call wtap2 (S20,b16nz,16,nz,4)
      call wtap2 (SN20,b5nz,5,nz,4)
      call wtap2 (S21,b16nz,16,nz,4)
      call wtap2 (SN21,b5nz,5,nz,4)
      call wtap2 (SCNT,b16nz,16,nz,4)
      call wtap2 (SNCNT,b5nz,5,nz,4)
      call wtap2 (SS9,b16nz,16,nz,4)
      call wtap2 (SSN9,b5nz,5,nz,4)
      call wtap2 (SS10,b16nz,16,nz,4)
      call wtap2 (SSN10,b5nz,5,nz,4)
      call wtap2 (SS11,b16nz,16,nz,4)
      call wtap2 (SSN11,b5nz,5,nz,4)
      call wtap2 (SS12,b16nz,16,nz,4)
      call wtap2 (SSN12,b5nz,5,nz,4)
      call wtap2 (SS13,b16nz,16,nz,4)
      call wtap2 (SSN13,b5nz,5,nz,4)
      call wtap2 (SS14,b16nz,16,nz,4)
      call wtap2 (SSN14,b5nz,5,nz,4)
      call wtap2 (SS15,b16nz,16,nz,4)
      call wtap2 (SSN15,b5nz,5,nz,4)
      call wtap2 (SS16,b16nz,16,nz,4)
      call wtap2 (SSN16,b5nz,5,nz,4)
      call wtap2 (SS17,b16nz,16,nz,4)
      call wtap2 (SSN17,b5nz,5,nz,4)
      call wtap2 (SS18,b16nz,16,nz,4)
      call wtap2 (SSN18,b5nz,5,nz,4)
      call wtap2 (SS19,b16nz,16,nz,4)
      call wtap2 (SSN19,b5nz,5,nz,4)
      call wtap2 (SS20,b16nz,16,nz,4)
      call wtap2 (SSN20,b5nz,5,nz,4)
      call wtap2 (SS21,b16nz,16,nz,4)
      call wtap2 (SSN21,b5nz,5,nz,4)
      call wtap2 (SSCNT,b16nz,16,nz,4)
      call wtap2 (SSNCNT,b5nz,5,nz,4)
      call wtap1 (tls,bnz,nz,4)
      call wtap1 (qls,bnz,nz,4)
      call wtap1 (cld,bnz,nz,4)
      call wtap1 (cldu,bnz,nz,4)
      call wtap1 (cldd,bnz,nz,4)
      call wtap1 (smf,bnz,nz,4)
      call wtap1 (smu,bnz,nz,4)
      call wtap1 (smd,bnz,nz,4)
      call wtap1 (smn,bnz,nz,4)
      call wtap1 (sft,bnz,nz,4)
      call wtap1 (sfq,bnz,nz,4)
      call wtap1 (stf1,bnz,nz,4)
      call wtap1 (stu1,bnz,nz,4)
      call wtap1 (std1,bnz,nz,4)
      call wtap1 (stn1,bnz,nz,4)
      call wtap1 (stf2,bnz,nz,4)
      call wtap1 (stu2,bnz,nz,4)
      call wtap1 (std2,bnz,nz,4)
      call wtap1 (stn2,bnz,nz,4)
      call wtap1 (stf,bnz,nz,4)
      call wtap1 (stu,bnz,nz,4)
      call wtap1 (std,bnz,nz,4)
      call wtap1 (stn,bnz,nz,4)
      call wtap1 (sqf1,bnz,nz,4)
      call wtap1 (squ1,bnz,nz,4)
      call wtap1 (sqd1,bnz,nz,4)
      call wtap1 (sqn1,bnz,nz,4)
      call wtap1 (sqf2,bnz,nz,4)
      call wtap1 (squ2,bnz,nz,4)
      call wtap1 (sqd2,bnz,nz,4)
      call wtap1 (sqn2,bnz,nz,4)
      call wtap1 (sqf,bnz,nz,4)
      call wtap1 (squ,bnz,nz,4)
      call wtap1 (sqd,bnz,nz,4)
      call wtap1 (sqn,bnz,nz,4)
      call wtap1 (sc,bnz,nz,4)
      call wtap1 (se,bnz,nz,4)
      call wtap1 (scu1,bnz,nz,4)
      call wtap1 (sed1,bnz,nz,4)
      call wtap1 (aub,bnz,nz,4)
      call wtap1 (avb,bnz,nz,4)
      call wtap1 (sut1,bnz,nz,4)
      call wtap1 (svt1,bnz,nz,4)
      call wtap1 (suu1,bnz,nz,4)
      call wtap1 (sud1,bnz,nz,4)
      call wtap1 (svu1,bnz,nz,4)
      call wtap1 (svd1,bnz,nz,4)
      call wtap1 (suu2,bnz,nz,4)
      call wtap1 (sud2,bnz,nz,4)
      call wtap1 (svu2,bnz,nz,4)
      call wtap1 (svd2,bnz,nz,4)
      call wtap1 (auu1,bnz,nz,4)
      call wtap1 (aud1,bnz,nz,4)
      call wtap1 (avu1,bnz,nz,4)
      call wtap1 (avd1,bnz,nz,4)
      call wtap1 (auu2,bnz,nz,4)
      call wtap1 (aud2,bnz,nz,4)
      call wtap1 (avu2,bnz,nz,4)
      call wtap1 (avd2,bnz,nz,4)
      call wtap3 (thom,bnz47,nz,4,7,4)
      call wtap3 (tdw,bnz47,nz,4,7,4)
      call wtap3 (tmlt,bnz47,nz,4,7,4)
      call wtap3 (saut,bnz47,nz,4,7,4)
      call wtap3 (saci,bnz47,nz,4,7,4)
      call wtap3 (sacw,bnz47,nz,4,7,4)
      call wtap3 (raci,bnz47,nz,4,7,4)
      call wtap3 (tacr,bnz47,nz,4,7,4)
      call wtap3 (raut,bnz47,nz,4,7,4)
      call wtap3 (racw,bnz47,nz,4,7,4)
      call wtap3 (sfw,bnz47,nz,4,7,4)
      call wtap3 (sfi,bnz47,nz,4,7,4)
      call wtap3 (gacs,bnz47,nz,4,7,4)
      call wtap3 (gacw,bnz47,nz,4,7,4)
      call wtap3 (gaci,bnz47,nz,4,7,4)
      call wtap3 (gacr,bnz47,nz,4,7,4)
      call wtap3 (gwet,bnz47,nz,4,7,4)
      call wtap3 (gaut,bnz47,nz,4,7,4)
      call wtap3 (racs,bnz47,nz,4,7,4)
      call wtap3 (sacr,bnz47,nz,4,7,4)
      call wtap3 (gfr,bnz47,nz,4,7,4)
      call wtap3 (smlt,bnz47,nz,4,7,4)
      call wtap3 (gmlt,bnz47,nz,4,7,4)
      call wtap3 (sdep,bnz47,nz,4,7,4)
      call wtap3 (ssub,bnz47,nz,4,7,4)
      call wtap3 (gsub,bnz47,nz,4,7,4)
      call wtap3 (pern,bnz47,nz,4,7,4)
      call wtap3 (d3ri,bnz47,nz,4,7,4)
      call wtap3 (d3ir,bnz47,nz,4,7,4)
      call wtap3 (d2sr,bnz47,nz,4,7,4)
      call wtap3 (d2rs,bnz47,nz,4,7,4)
      call wtap3 (gdry,bnz47,nz,4,7,4)
      call wtap3 (coc,bnz47,nz,4,7,4)
      call wtap3 (coe,bnz47,nz,4,7,4)
      call wtap3 (smf0,bnz47,nz,4,7,4)
      call wtap3 (qc0,bnz47,nz,4,7,4)
      call wtap3 (qr0,bnz47,nz,4,7,4)
      call wtap3 (qi0,bnz47,nz,4,7,4)
      call wtap3 (qs0,bnz47,nz,4,7,4)
      call wtap3 (qg0,bnz47,nz,4,7,4)
      call wtap3 (sqc0,bnz47,nz,4,7,4)
      call wtap3 (sqr0,bnz47,nz,4,7,4)
      call wtap3 (sqi0,bnz47,nz,4,7,4)
      call wtap3 (sqs0,bnz47,nz,4,7,4)
      call wtap3 (sqg0,bnz47,nz,4,7,4)
      call wtap3 (erns,bnz47,nz,4,7,4)
      call wtap3 (wgrs,bnz47,nz,4,7,4)
      call wtap3 (qsws,bnz47,nz,4,7,4)
      call wtap3 (srsw,bnz47,nz,4,7,4)
      call wtap3 (srlw,bnz47,nz,4,7,4)
      call wtap3 (tut1,bnz47,nz,4,7,4)
      call wtap3 (tut2,bnz47,nz,4,7,4)
      call wtap3 (tvt1,bnz47,nz,4,7,4)
      call wtap3 (tvt2,bnz47,nz,4,7,4)
      call wtap3 (tstf,bnz47,nz,4,7,4)
      call wtap3 (tstf1,bnz47,nz,4,7,4)
      call wtap3 (tstf2,bnz47,nz,4,7,4)
      call wtap3 (tsqf,bnz47,nz,4,7,4)
      call wtap3 (tsqf1,bnz47,nz,4,7,4)
      call wtap3 (tsqf2,bnz47,nz,4,7,4)
      call wtap3 (tsqq,bnz47,nz,4,7,4)
      call wtap3 (tsqq1,bnz47,nz,4,7,4)
      call wtap3 (fcld,bnz47,nz,4,7,4)
      CALL WTAP3 (OTHERT_ADD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (OTHERQ_ADD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNTH,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNTV,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNTD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNQH,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNQV,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNQD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNHH,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNHV,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNHD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (T_BM,BNZ47,NZ,4,7,4)
      CALL WTAP3 (Q_BM,BNZ47,NZ,4,7,4)
c      call wtap3 (ceds1i,bnxnz4,nx,nz,4,4)
c      call wtap3 (ceds2i,bnxnz4,nx,nz,4,4)
c      call wtap3 (tstfi,bnxnz4,nx,nz,4,4)
c      call wtap3 (tsqfi,bnxnz4,nx,nz,4,4)
c      call wtap3 (rsli,bnxnz4,nx,nz,4,4)
      call wtap2 (tb0,bnz4,nz,4,4)
      call wtap2 (qb0,bnz4,nz,4,4)
c      print *,'wtapnew'
c      call wtap1 (avett,bnnt,nnt,4)
c      call wtap1 (avetq,bnnt,nnt,4)
c      call wtap1 (avet950,bnnt,nnt,4)
c      call wtap1 (acoc,bnz,nz,4)
c      call wtap1 (acoe,bnz,nz,4)
c      call wtap1 (agaut,bnz,nz,4)
c      call wtap1 (asdep,bnz,nz,4)
c      call wtap1 (agsub,bnz,nz,4)
c      call wtap1 (apern,bnz,nz,4)
c      call wtap1 (aqls,bnz,nz,4)
c      call wtap1 (altrans,bnz,nz,4)
c      call wtap1 (actrans,bnz,nz,4)
c      call wtap1 (total_mq,bnz,nz,4)
c      call wtap1 (agwet,bnz,nz,4)
c      call wtap1 (afgr,bnz,nz,4)
c      call wtap1 (avelw,bnz,nz,4)
c      call wtap1 (avesw,bnz,nz,4)
c      call wtap1 (amelt,bnz,nz,4)
c      call wtap1 (aother,bnz,nz,4)
c      call wtap1 (total_mt,bnz,nz,4)
c      call wtap1 (tavet,bnz,nz,4)
c      call wtap1 (qavet,bnz,nz,4)
c      call wtap1 (t_sfcq,bnnt,nnt,4)
c      call wtap1 (t_sfct,bnnt,nnt,4)
c      call wtap1 (t_ub,bnnt,nnt,4)
c      call wtap1 (t_microq,bnnt,nnt,4)
c      call wtap1 (t_largeq,bnnt,nnt,4)
c      call wtap1 (t_larget,bnnt,nnt,4)
c      call wtap1 (t_microt,bnnt,nnt,4)
c      call wtap1 (tsfc_1,bnx,nx,4)
c      call wtap1 (qsfc_1,bnx,nx,4)
c      call wtap1 (tlsw,bnz,nz,4)
c      call wtap1 (qlsw,bnz,nz,4)
c      call wtap1 (ttlsw,bnz,nz,4)
c      call wtap1 (tqlsw,bnz,nz,4)
c      call wtap1 (tcof,bnz,nz,4)
c      call wtap1 (textra,bnz,nz,4)
c      call wtap1 (qextra,bnz,nz,4)
c      call wtap1 (tlsw1,bnz,nz,4)
c      call wtap1 (qlsw1,bnz,nz,4)
c      call wtap1 (ttlsw1,bnz,nz,4)
c      call wtap1 (tqlsw1,bnz,nz,4)

       CALL WTAP_UV_S

c     ***********************************
 2222 if(isec.EQ.isec/imp*imp) THEN
        write(6,201) n,aminut,smt,smv,smq,smc,sme,smr
        write(6,202) n,aminut,smqc,smqr,smqi,smqs,smqg
c
        write(6,207) n,aminut,smrc
c
 7207 format(10h time step,i5,6h time=,f8.2,9h min dep=,e11.4,5h sub=,
     1  e11.4,5h qrs=,e11.4,5h qrl=,e11.4,5h mlt=,e11.4,5h frz=,e11.4)
        write(6,7207) n,aminut,sddd,ssss,shhh,sccc,smmm,sfff
c
        write(6,1234) aco,aan,aco1,aan1
        write(6,2345) lconv,lanvl,lnspt
        ICOT=4
        IF (NXI .EQ. 128) ICOT=1
        IF (NXI .EQ. 256) ICOT=2
        IF (NXI .EQ. 1024) ICOT=8
        write(6,3456) (iv(i),i=2,iles,ICOT)
        do 33 k=2,kles
          do 33 i=2,iles
            aaa1(i,k)=blk
            y1(i)=.25*(qcl(i,k)+qcl(i,k-1)+qcl(i-1,k)+qcl(i-1,k-1))
            y2(i)=.25*(qrn(i,k)+qrn(i,k-1)+qrn(i-1,k)+qrn(i-1,k-1))
            if(y1(i).ge.cotor) aaa1(i,k)=dot
            if(y2(i).ge.cotor) aaa1(i,k)=exm
            if(y1(i).ge.cotor.and.y2(i).ge.cotor) aaa1(i,k)=com
   33   continue
          do 34 k1=2,kles
            k=kles+2-k1
            km=k-1
   34     write(6,203) km,(aaa1(i,k),i=2,iles,ICOT)
        IF (ICE .EQ. 1) THEN
          do 35 k=2,kles
            do 35 i=2,iles
              aaa1(i,k)=blk
              y1(i)=.25*(qci(i,k)+qci(i,k-1)+qci(i-1,k)+qci(i-1,k-1))
              y2(i)=.25*(qcs(i,k)+qcs(i,k-1)+qcs(i-1,k)+qcs(i-1,k-1))
              if (y1(i) .ge. cotor) aaa1(i,k)=pice
              if (y2(i) .ge. cotor) aaa1(i,k)=ago
              if (y1(i) .ge. cotor .and. y2(i) .ge. cotor) aaa1(i,k)=abt
   35     continue
          do 36 k1=2,kles
            k=kles+2-k1
            km=k-1
   36     write(6,203) km,(aaa1(i,k),i=2,iles,ICOT)
          do 37 k=2,kles
            do 37 i=2,iles
              aaa1(i,k)=blk
              y1(i)=.25*(qcs(i,k)+qcs(i,k-1)+qcs(i-1,k)+qcs(i-1,k-1))
              y2(i)=.25*(qcg(i,k)+qcg(i,k-1)+qcg(i-1,k)+qcg(i-1,k-1))
              if (y1(i) .ge. cotor) aaa1(i,k)=ast
              if (y2(i) .ge. cotor) aaa1(i,k)=ago
              if (y1(i) .ge. cotor .and. y2(i) .ge. cotor) aaa1(i,k)=abt
   37     continue
          do 38 k1=2,kles
            k=kles+2-k1
            km=k-1
   38     write(6,203) km,(aaa1(i,k),i=2,iles,ICOT)
        ENDIF
      ENDIF
c     ***********************************
      if (isec.eq.isec/120*120) then
         call tmax (w,'W   ',aminut,wmax,wmin,1)
      endif
      if (isec.eq.isec/imx*imx) then
        m=m+1
        lc=0
        ls=0
        ln=0
        do 21 i=2,iles
          lc=lc+ics(i,2)
          ls=ls+ics(i,3)
   21     ln=ln+ics(i,4)
        cnv(m)=float(lc)
        anv(m)=float(ls)
        snp(m)=float(ln)
        call tmax (w,'W   ',aminut,wmax,wmin,1)
        rwmax(m)=wmax
        rwmin(m)=wmin
      endif
      if (isec.eq.isec/imx1*imx1) then
         call tmax (dpt,'DPT ',aminut,wmax,wmin,2)
         call tmax (dqv,'DQV ',aminut,wmax,wmin,3)
         call tmax (qcl,'QCL ',aminut,wmax,wmin,3)
         call tmax (qrn,'QRN ',aminut,wmax,wmin,3)
         IF (ICE .EQ. 1) THEN
            call tmax (qci,'QCI ',aminut,wmax,wmin,3)
            call tmax (qcs,'QCS ',aminut,wmax,wmin,3)
            call tmax (qcg,'QCG ',aminut,wmax,wmin,3)
         ENDIF
         call tmax (ak,'AK  ',aminut,wmax,wmin,4)
       endif
       if (isec.eq.isec/imx*imx) then
         m1=m1+1
         tqc(m1)=smqc
         tqr(m1)=smqr
         tqi(m1)=smqi
         tqs(m1)=smqs
         tqg(m1)=smqg
c
         call ptpb(m1)
C
          DO I=2,ILES
            IF (IRICONT1 .EQ. 0) THEN
                RI120AR(I)=0.
            ELSE
                RI120AR(I)=RI120AR(I)/FLOAT(IRICONT1)
            ENDIF
          ENDDO
         DO 39 I=1,NXI
            IH=I+1
           RSVAR(M1,I)=RI120AR(IH)
           RST(M1,I)=RI(IH)
   39      RSV(M1,I)=IV(IH)
          IRICONT1=0
           DO I=2,ILES
             RI120AR(I)=0.
           ENDDO
       endif
c     *********************************
      if(isec.eq.isec/21600*21600) then
      write(6,125)
      do 81 k1=2,kles
       k=kles+2-k1
       km=k-1
   81  write(6,134) km,st(k),sv(k),sq(k),sc(k),se(k),sqa(k),cld(k),
     1    cldu(k),cldd(k),ub(k),vb(k),wb(k)
c
      write(6,5125)
      do k1=2,kles
         k=kles+2-k1
       write(6,134) k,s_dep(k),s_sub(k),s_qrs(k),s_qrl(k),s_mel(k),
     1              s_frz(k)
      enddo
 5125 FORMAT(1X,2HLV,6X,3HDEP,6X,3HSUB,6X,3HQRH,6X,3HQRC,6X,3HMET,6X,3HF
     1EZ)
c
      write(6,1025)
      do k1=2,kles
        k=kles+2-k1
        km=k-1
       write(6,134) km,Q1ZT(k),Q2ZT(k),Q1T(k),Q2T(k),aq1t(k),aq2t(k)
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
13666 format(111h lv       ut       vt       wt      q1t      q2t      q
     11z       q2z                                            )

      write(6,13666)

       do k=2,kles
         y1(k)=0.
         y2(k)=0.
         y3(k)=0.
         y4(k)=0.
         y5(k)=0.
         y6(k)=0.
       enddo
       do k=2,kles
       do i=2,iles
         y1(k)=y1(k)+ww1(i,k)
         y2(k)=y2(k)+uu1(i,k)
         y3(k)=y3(k)+vv1(i,k)
         y4(k)=y4(k)+w(i,k)
         y5(k)=y5(k)+u(i,k)
         y6(k)=y6(k)+v(i,k)
       enddo
       enddo
      do k1=2,kles
       k=kles+2-k1
       km=k-1
        a1=ubt(k)*6.*3600.
        a2=vbt(k)*6.*3600.
        a3=wbt(k)*6.*3600.
        a4=q1t(k)*6.*3600.
        a5=q2t(k)*6.*3600.
        a6=q1zt(k)*6.*3600.
        a7=q2zt(k)*6.*3600.
       write(6,134) km,a1,a2,a3,a4,a5,a6,a7,y1(k)*ril2,y2(k)*ril2,
     1              y3(k)*ril2,y4(k)*ril2,y5(k)*ril2,y6(k)*ril2
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(6,126)
      do 82 k1=2,kles
       k=kles+2-k1
       km=k-1
   82  write(6,134) km,sqc(k),sqr(k),sqi(k),sqs(k),sqg(k),stqc(k),
     1    stqr(k),stqi(k),stqs(k),stqg(k)
      do 83 k=2,kles
       stv(k)=stf(k)+.61*tb(k)*sqf(k)
       sft(k)=am(k)*(stf(k)-stf(k+1))/(rho(k)*dz)
   83  sfq(k)=am(k)*(sqf(k)-sqf(k+1))/(rho(k)*dz)
      write(6,139)
      do 84 k1=2,kles
       k=kles+2-k1
   84  write(6,134) k,smf(k),smu(k),smd(k),stf(k),stu(k),std(k),
     1   sqf(k),squ(k),sqd(k),sdt(k),sdq(k),stv(k),sft(k),sfq(k)
      endif
c     ********************************
 5555 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if (isec .eq. isec/itime_ave*itime_ave) then
c         mmave=mmave+1
c         call dom_ave (aminut,tave,denq,thee)          
ccccccccc
c         dt_stat=ndt_stat
c         call budget_mt (dt_stat,total_mqt,total_mtt)
c         call var_ave (total_mqt,tqqq2)
c         call var_ave (qls,tqqq3)
c         call var_avet (total_mtt,tttt2)
c         call var_avet (tls,tttt3)
ccccccccc
c         avett(mmave)=tave
c         avetq(mmave)=denq
c         avet950(mmave)=thee
c         t_sfcq(mmave)=sqfff2
c         t_sfct(mmave)=stfff2/denave
c         t_ub(mmave)=ub1(2)
c         t_microq(mmave)=tqqq2
c         t_largeq(mmave)=tqqq3
c         t_larget(mmave)=tttt3
c         t_microt(mmave)=tttt2
c      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ISMG=6
      IADVH=IADVH1
        IF (IJKADV .EQ. 1) IADVH=2
        irf=1
      CALL AKCOEF (PI)
c      if (ijkadv .eq. 1) then
        DO 900 K=1,KMAX
        DO 900 I=1,IMAX
          UMD(I,K)=UU1(I,K)
          vmd(i,k)=vv1(i,k)
  900     WMD(I,K)=WW1(I,K)
c      endif
       irf=1
       IADVH=IADVH1
c
       call RINIT_UV2
c
      call advwuv (nudge)
      call cmpf
      call slvpi (1)
        id_uv=0
c         if (isec .eq. isec/120*120) id_uv=1
          if (isec .eq. isec/900*900) id_uv=1 ! 8/13/01
      call cmpwuv (id_uv)
c     ***********************************
      if(isec.eq.isec/240*240) then
        ARF=XXW(2,KLes)
        IF (LIPPS .EQ. 1) ARF=XXW(2,KLes)/TB(KLes)
        DO 502 K=2,KLES
        DO 502 I=2,ILES
          IF (LIPPS .EQ. 1) THEN
            AAA(I,K)=XXW(I,K)/TB(K)-ARF
          ELSE
            AAA(I,K)=XXW(I,K)-ARF
          ENDIF
  502   CONTINUE
        call tmax (aaa,'PI  ',aminut,wmax,wmin,3)
       endif
c     ***********************************
      if(isec.lt.iwtp) go to 9999
      if(isec.ne.isec/iwp*iwp) go to 9999
      ARF=XXW(2,KLes)
      IF (LIPPS .EQ. 1) ARF=XXW(2,KLes)/TB(KLes)
      DO 5052 K=2,KLES
      DO 5052 I=2,ILES
       IF (LIPPS .EQ. 1) THEN
         AAA(I,K)=XXW(I,K)/TB(K)-ARF
       ELSE
         AAA(I,K)=XXW(I,K)-ARF
       ENDIF
 5052 CONTINUE
        call wtap2 (aaa,bnxnz,nx,nz,4)
c     ************************
 9999 if(isec.lt.iend) go to 3333
c     ************************
      do 80 l=4,4
      do 80 j=1,4
      write(6,133) l,j,lc,ls,ln
      write(6,135)
      do 85 k1=2,kles
       k=kles+2-k1
   85  write(6,134) k,thom(k,j,l),tdw(k,j,l),tmlt(k,j,l),saut(k,j,l),
     1   saci(k,j,l),sacw(k,j,l),raci(k,j,l),tacr(k,j,l),raut(k,j,l),
     2   racw(k,j,l),sfw(k,j,l),sfi(k,j,l),gacs(k,j,l),gacw(k,j,l)
      write(6,136)
      do 86 k1=2,kles
       k=kles+2-k1
   86  write(6,134) k,gaci(k,j,l),gacr(k,j,l),gwet(k,j,l),gaut(k,j,l),
     1   racs(k,j,l),sacr(k,j,l),gfr(k,j,l),smlt(k,j,l),gmlt(k,j,l),
     2   sdep(k,j,l),ssub(k,j,l),gsub(k,j,l)
      write(6,137)
      do 87 k1=2,kles
       k=kles+2-k1
   87  write(6,134) k,pern(k,j,l),d3ri(k,j,l),d3ir(k,j,l),d2sr(k,j,l),
     1   d2rs(k,j,l),gdry(k,j,l),erns(k,j,l),wgrs(k,j,l),qsws(k,j,l),
     2   srsw(k,j,l),srlw(k,j,l)
      write(6,138)
      do 89 k1=2,kles
       k=kles+2-k1
   89  write(6,134) k,coc(k,j,l),coe(k,j,l),smf0(k,j,l),qc0(k,j,l),
     1   qr0(k,j,l),qi0(k,j,l),qs0(k,j,l),qg0(k,j,l),sqc0(k,j,l),
     2   sqr0(k,j,l),sqi0(k,j,l),sqs0(k,j,l),sqg0(k,j,l)
      write(6,140)
      do 96 k1=2,kles
       k=kles+2-k1
   96  write(6,134) k,tut1(k,j,l),tut2(k,j,l),tvt1(k,j,l),tvt2(k,j,l),
     1   tstf(k,j,l),tstf1(k,j,l),tstf2(k,j,l),tsqf(k,j,l),tsqf1(k,j,l),
     2   tsqf2(k,j,l),tsqq(k,j,l),tsqq1(k,j,l)
   80 continue
cccccccccccccccccccccccccccc
       write(6,1234) aco,aan
      write(6,160)
      do 991 i=1,nxi
      do 91 ntime=1,80
   91  iv(ntime)=int (rst(ntime,i))
  991  write(6,204) i,(iv(m),m=1,80,2)
      write(6,160)
      do 993 i=1,nxi
      do 93 ntime=1,80
   93  jv(ntime)=int (rsv(ntime,i))
  993  write(6,204) i,(jv(m),m=1,80,2)
      write(6,160)
      do 992 i=1,nxi
      do 92 ntime=81,nt
   92  iv(ntime)=int (rst(ntime,i))
  992  write(6,204) i,(iv(m),m=81,120)
      write(6,160)
      do 994 i=1,nxi
      do 94 ntime=81,nt
   94  jv(ntime)=int (rsv(ntime,i))
  994  write(6,204) i,(jv(m),m=81,120)
       call RAINSTAT (RST)
c    *************************************
      call wtap2 (rst,bntnxi,nt,nxi,4)
      call wtap2 (pcltop,bntnxi,nt,nxi,4)
      call wtap2 (pclbot,bntnxi,nt,nxi,4)
      call wtap2 (rsv,bntnxi,nt,nxi,4)
      call wtap2 (rsvar,bntnxi,nt,nxi,4)
      call wtap1 (tqc,bnt,nt,4)
      call wtap1 (tqr,bnt,nt,4)
      call wtap1 (tqi,bnt,nt,4)
      call wtap1 (tqs,bnt,nt,4)
      call wtap1 (tqg,bnt,nt,4)
      call wtap1 (rwmax,bnt3,nt3,4)
      call wtap1 (rwmin,bnt3,nt3,4)
      call wtap1 (cnv,bnt3,nt3,4)
      call wtap1 (anv,bnt3,nt3,4)
      call wtap1 (snp,bnt3,nt3,4)
   
      call wtap2 (ak,bnxnz,nx,nz,4)
      call wtap2 (dpt,bnxnz,nx,nz,4)
      call wtap2 (dqv,bnxnz,nx,nz,4)
      call wtap2 (qcl,bnxnz,nx,nz,4)
      call wtap2 (qrn,bnxnz,nx,nz,4)
      call wtap2 (qci,bnxnz,nx,nz,4)
      call wtap2 (qcs,bnxnz,nx,nz,4)
      call wtap2 (qcg,bnxnz,nx,nz,4)
      call wtap2 (u,bnxnz,nx,nz,4)
      call wtap2 (v,bnxnz,nx,nz,4)
      call wtap2 (w,bnxnz,nx,nz,4)
      call wtap2 (rsw,bnxnz,nx,nz,4)
      call wtap2 (rlw,bnxnz,nx,nz,4)

       CALL WTAP_UV
       CALL WTAP_UV_S

      if(icps_bm .eq. 1) then
      CALL WTAP2 (TMODO,BNXNZ,NX,NZ,4)
      CALL WTAP2 (QMODO,BNXNZ,NX,NZ,4)
      endif
      call wtap1 (ub,bnz,nz,4)
      call wtap1 (vb,bnz,nz,4)
      call wtap1 (ub1,bnz,nz,4)
      call wtap1 (vb1,bnz,nz,4)
      call wtap1 (tb,bnz,nz,4)
      call wtap1 (qb,bnz,nz,4)
      call wtap1 (rho,bnz,nz,4)
      call wtap1 (rho1,bnz,nz,4)
      call wtap1 (ta,bnz,nz,4)
      call wtap1 (qa,bnz,nz,4)
      call wtap1 (ta1,bnz,nz,4)
      call wtap1 (qa1,bnz,nz,4)
      call wtap1 (pi,bnz,nz,4)
      call wtap1 (p0,bnz,nz,4)
      call wtap1 (wb,bnz,nz,4)
      call wtap1 (wbt,bnz,nz,4)
      call wtap1 (fd,bnz,nz,4)
      call wtap1 (fe,bnz,nz,4)
      call wtap1 (sqc,bnz,nz,4)
      call wtap1 (sqr,bnz,nz,4)
      call wtap1 (sqi,bnz,nz,4)
      call wtap1 (sqs,bnz,nz,4)
      call wtap1 (sqg,bnz,nz,4)
      call wtap1 (stqc,bnz,nz,4)
      call wtap1 (stqr,bnz,nz,4)
      call wtap1 (stqi,bnz,nz,4)
      call wtap1 (stqs,bnz,nz,4)
      call wtap1 (stqg,bnz,nz,4)
      call wtap1 (aq1t,bnz,nz,4)
      call wtap1 (aq2t,bnz,nz,4)
      call wtap1 (aq1zt,bnz,nz,4)
      call wtap1 (aq2zt,bnz,nz,4)
      call wtap1 (suw,bnx,nx,4)
      call wtap1 (svw,bnx,nx,4)
      call wtap1 (swt,bnx,nx,4)
      call wtap1 (swq,bnx,nx,4)
      call wtap1 (tairsfc,bnx,nx,4)
      call wtap1 (qairsfc,bnx,nx,4)
      call wtap1 (thairsf,bnx,nx,4)
      call wtap1 (pairsfc,bnx,nx,4)
      call wtap2 (sun_4,bnx4,nx,4,4)
      call wtap2 (rflux,bnx7,nx,7,4)

c      call wtap2(q1a_g_h,bnxnz,nx,nz,4)
c      call wtap2(q1a_g_v,bnxnz,nx,nz,4)
c      call wtap2(q1a_d_h,bnxnz,nx,nz,4)
c      call wtap2(q1a_d_v,bnxnz,nx,nz,4)
c      call wtap2(q2a_g_h,bnxnz,nx,nz,4)
c      call wtap2(q2a_g_v,bnxnz,nx,nz,4)
c      call wtap2(q2a_d_h,bnxnz,nx,nz,4)
c      call wtap2(q2a_d_v,bnxnz,nx,nz,4)
c      call wtap2(q1a_hyd,bnxnz,nx,nz,4)
c      call wtap2(q2a_hyd,bnxnz,nx,nz,4)
c      call wtap2(q1a_rad,bnxnz,nx,nz,4)

      call wtap0 (smr,4)
      call wtap0 (smrc,4)
      call wtap1 (ri,bnx,nx,4)
      call wtap1 (ar_c,bnx,nx,4)
      call wtap1 (ar,bnx,nx,4)
      call wtap1 (ri180,bnx,nx,4)
      call wtap1 (ri120ar,bnx,nx,4)

      if(istaton.eq.1)then

      call wtap3 (S1,b16nz14,16,nz,14,4)
      call wtap3 (SN1,b5nz14,5,nz,14,4)
      call wtap2 (S9,b16nz,16,nz,4)
      call wtap2 (SN9,b5nz,5,nz,4)
      call wtap2 (S10,b16nz,16,nz,4)
      call wtap2 (SN10,b5nz,5,nz,4)
      call wtap2 (S11,b16nz,16,nz,4)
      call wtap2 (SN11,b5nz,5,nz,4)
      call wtap2 (S12,b16nz,16,nz,4)
      call wtap2 (SN12,b5nz,5,nz,4)
      call wtap2 (S13,b16nz,16,nz,4)
      call wtap2 (SN13,b5nz,5,nz,4)
      call wtap2 (S14,b16nz,16,nz,4)
      call wtap2 (SN14,b5nz,5,nz,4)
      call wtap2 (S15,b16nz,16,nz,4)
      call wtap2 (SN15,b5nz,5,nz,4)
      call wtap2 (S16,b16nz,16,nz,4)
      call wtap2 (SN16,b5nz,5,nz,4)
      call wtap2 (S17,b16nz,16,nz,4)
      call wtap2 (SN17,b5nz,5,nz,4)
      call wtap2 (S18,b16nz,16,nz,4)
      call wtap2 (SN18,b5nz,5,nz,4)
      call wtap2 (S19,b16nz,16,nz,4)
      call wtap2 (SN19,b5nz,5,nz,4)
      call wtap2 (S20,b16nz,16,nz,4)
      call wtap2 (SN20,b5nz,5,nz,4)
      call wtap2 (S21,b16nz,16,nz,4)
      call wtap2 (SN21,b5nz,5,nz,4)
      call wtap2 (SCNT,b16nz,16,nz,4)
      call wtap2 (SNCNT,b5nz,5,nz,4)
      call wtap2 (SS9,b16nz,16,nz,4)
      call wtap2 (SSN9,b5nz,5,nz,4)
      call wtap2 (SS10,b16nz,16,nz,4)
      call wtap2 (SSN10,b5nz,5,nz,4)
      call wtap2 (SS11,b16nz,16,nz,4)
      call wtap2 (SSN11,b5nz,5,nz,4)
      call wtap2 (SS12,b16nz,16,nz,4)
      call wtap2 (SSN12,b5nz,5,nz,4)
      call wtap2 (SS13,b16nz,16,nz,4)
      call wtap2 (SSN13,b5nz,5,nz,4)
      call wtap2 (SS14,b16nz,16,nz,4)
      call wtap2 (SSN14,b5nz,5,nz,4)
      call wtap2 (SS15,b16nz,16,nz,4)
      call wtap2 (SSN15,b5nz,5,nz,4)
      call wtap2 (SS16,b16nz,16,nz,4)
      call wtap2 (SSN16,b5nz,5,nz,4)
      call wtap2 (SS17,b16nz,16,nz,4)
      call wtap2 (SSN17,b5nz,5,nz,4)
      call wtap2 (SS18,b16nz,16,nz,4)
      call wtap2 (SSN18,b5nz,5,nz,4)
      call wtap2 (SS19,b16nz,16,nz,4)
      call wtap2 (SSN19,b5nz,5,nz,4)
      call wtap2 (SS20,b16nz,16,nz,4)
      call wtap2 (SSN20,b5nz,5,nz,4)
      call wtap2 (SS21,b16nz,16,nz,4)
      call wtap2 (SSN21,b5nz,5,nz,4)
      call wtap2 (SSCNT,b16nz,16,nz,4)
      call wtap2 (SSNCNT,b5nz,5,nz,4)
      call wtap1 (tls,bnz,nz,4)
      call wtap1 (qls,bnz,nz,4)
      call wtap1 (cld,bnz,nz,4)
      call wtap1 (cldu,bnz,nz,4)
      call wtap1 (cldd,bnz,nz,4)
      call wtap1 (smf,bnz,nz,4)
      call wtap1 (smu,bnz,nz,4)
      call wtap1 (smd,bnz,nz,4)
      call wtap1 (smn,bnz,nz,4)
      call wtap1 (sft,bnz,nz,4)
      call wtap1 (sfq,bnz,nz,4)
      call wtap1 (stf1,bnz,nz,4)
      call wtap1 (stu1,bnz,nz,4)
      call wtap1 (std1,bnz,nz,4)
      call wtap1 (stn1,bnz,nz,4)
      call wtap1 (stf2,bnz,nz,4)
      call wtap1 (stu2,bnz,nz,4)
      call wtap1 (std2,bnz,nz,4)
      call wtap1 (stn2,bnz,nz,4)
      call wtap1 (stf,bnz,nz,4)
      call wtap1 (stu,bnz,nz,4)
      call wtap1 (std,bnz,nz,4)
      call wtap1 (stn,bnz,nz,4)
      call wtap1 (sqf1,bnz,nz,4)
      call wtap1 (squ1,bnz,nz,4)
      call wtap1 (sqd1,bnz,nz,4)
      call wtap1 (sqn1,bnz,nz,4)
      call wtap1 (sqf2,bnz,nz,4)
      call wtap1 (squ2,bnz,nz,4)
      call wtap1 (sqd2,bnz,nz,4)
      call wtap1 (sqn2,bnz,nz,4)
      call wtap1 (sqf,bnz,nz,4)
      call wtap1 (squ,bnz,nz,4)
      call wtap1 (sqd,bnz,nz,4)
      call wtap1 (sqn,bnz,nz,4)
      call wtap1 (sc,bnz,nz,4)
      call wtap1 (se,bnz,nz,4)
      call wtap1 (scu1,bnz,nz,4)
      call wtap1 (sed1,bnz,nz,4)
      call wtap1 (aub,bnz,nz,4)
      call wtap1 (avb,bnz,nz,4)
      call wtap1 (sut1,bnz,nz,4)
      call wtap1 (svt1,bnz,nz,4)
      call wtap1 (suu1,bnz,nz,4)
      call wtap1 (sud1,bnz,nz,4)
      call wtap1 (svu1,bnz,nz,4)
      call wtap1 (svd1,bnz,nz,4)
      call wtap1 (suu2,bnz,nz,4)
      call wtap1 (sud2,bnz,nz,4)
      call wtap1 (svu2,bnz,nz,4)
      call wtap1 (svd2,bnz,nz,4)
      call wtap1 (auu1,bnz,nz,4)
      call wtap1 (aud1,bnz,nz,4)
      call wtap1 (avu1,bnz,nz,4)
      call wtap1 (avd1,bnz,nz,4)
      call wtap1 (auu2,bnz,nz,4)
      call wtap1 (aud2,bnz,nz,4)
      call wtap1 (avu2,bnz,nz,4)
      call wtap1 (avd2,bnz,nz,4)
      call wtap3 (thom,bnz47,nz,4,7,4)
      call wtap3 (tdw,bnz47,nz,4,7,4)
      call wtap3 (tmlt,bnz47,nz,4,7,4)
      call wtap3 (saut,bnz47,nz,4,7,4)
      call wtap3 (saci,bnz47,nz,4,7,4)
      call wtap3 (sacw,bnz47,nz,4,7,4)
      call wtap3 (raci,bnz47,nz,4,7,4)
      call wtap3 (tacr,bnz47,nz,4,7,4)
      call wtap3 (raut,bnz47,nz,4,7,4)
      call wtap3 (racw,bnz47,nz,4,7,4)
      call wtap3 (sfw,bnz47,nz,4,7,4)
      call wtap3 (sfi,bnz47,nz,4,7,4)
      call wtap3 (gacs,bnz47,nz,4,7,4)
      call wtap3 (gacw,bnz47,nz,4,7,4)
      call wtap3 (gaci,bnz47,nz,4,7,4)
      call wtap3 (gacr,bnz47,nz,4,7,4)
      call wtap3 (gwet,bnz47,nz,4,7,4)
      call wtap3 (gaut,bnz47,nz,4,7,4)
      call wtap3 (racs,bnz47,nz,4,7,4)
      call wtap3 (sacr,bnz47,nz,4,7,4)
      call wtap3 (gfr,bnz47,nz,4,7,4)
      call wtap3 (smlt,bnz47,nz,4,7,4)
      call wtap3 (gmlt,bnz47,nz,4,7,4)
      call wtap3 (sdep,bnz47,nz,4,7,4)
      call wtap3 (ssub,bnz47,nz,4,7,4)
      call wtap3 (gsub,bnz47,nz,4,7,4)
      call wtap3 (pern,bnz47,nz,4,7,4)
      call wtap3 (d3ri,bnz47,nz,4,7,4)
      call wtap3 (d3ir,bnz47,nz,4,7,4)
      call wtap3 (d2sr,bnz47,nz,4,7,4)
      call wtap3 (d2rs,bnz47,nz,4,7,4)
      call wtap3 (gdry,bnz47,nz,4,7,4)
      call wtap3 (coc,bnz47,nz,4,7,4)
      call wtap3 (coe,bnz47,nz,4,7,4)
      call wtap3 (smf0,bnz47,nz,4,7,4)
      call wtap3 (qc0,bnz47,nz,4,7,4)
      call wtap3 (qr0,bnz47,nz,4,7,4)
      call wtap3 (qi0,bnz47,nz,4,7,4)
      call wtap3 (qs0,bnz47,nz,4,7,4)
      call wtap3 (qg0,bnz47,nz,4,7,4)
      call wtap3 (sqc0,bnz47,nz,4,7,4)
      call wtap3 (sqr0,bnz47,nz,4,7,4)
      call wtap3 (sqi0,bnz47,nz,4,7,4)
      call wtap3 (sqs0,bnz47,nz,4,7,4)
      call wtap3 (sqg0,bnz47,nz,4,7,4)
      call wtap3 (erns,bnz47,nz,4,7,4)
      call wtap3 (wgrs,bnz47,nz,4,7,4)
      call wtap3 (qsws,bnz47,nz,4,7,4)
      call wtap3 (srsw,bnz47,nz,4,7,4)
      call wtap3 (srlw,bnz47,nz,4,7,4)
      call wtap3 (tut1,bnz47,nz,4,7,4)
      call wtap3 (tut2,bnz47,nz,4,7,4)
      call wtap3 (tvt1,bnz47,nz,4,7,4)
      call wtap3 (tvt2,bnz47,nz,4,7,4)
      call wtap3 (tstf,bnz47,nz,4,7,4)
      call wtap3 (tstf1,bnz47,nz,4,7,4)
      call wtap3 (tstf2,bnz47,nz,4,7,4)
      call wtap3 (tsqf,bnz47,nz,4,7,4)
      call wtap3 (tsqf1,bnz47,nz,4,7,4)
      call wtap3 (tsqf2,bnz47,nz,4,7,4)
      call wtap3 (tsqq,bnz47,nz,4,7,4)
      call wtap3 (tsqq1,bnz47,nz,4,7,4)
      call wtap3 (fcld,bnz47,nz,4,7,4)
      CALL WTAP3 (OTHERT_ADD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (OTHERQ_ADD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNTH,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNTV,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNTD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNQH,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNQV,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNQD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNHH,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNHV,BNZ47,NZ,4,7,4)
      CALL WTAP3 (SNHD,BNZ47,NZ,4,7,4)
      CALL WTAP3 (T_BM,BNZ47,NZ,4,7,4)
      CALL WTAP3 (Q_BM,BNZ47,NZ,4,7,4)
c      call wtap3 (ceds1i,bnxnz4,nx,nz,4,4)
c      call wtap3 (ceds2i,bnxnz4,nx,nz,4,4)
c      call wtap3 (tstfi,bnxnz4,nx,nz,4,4)
c      call wtap3 (tsqfi,bnxnz4,nx,nz,4,4)
c      call wtap3 (rsli,bnxnz4,nx,nz,4,4)
      call wtap2 (tb0,bnz4,nz,4,4)
      call wtap2 (qb0,bnz4,nz,4,4)
c      call wtap1 (avett,bnnt,nnt,4)
c      call wtap1 (avetq,bnnt,nnt,4)
c      call wtap1 (avet950,bnnt,nnt,4)
c      call wtap1 (acoc,bnz,nz,4)
c      call wtap1 (acoe,bnz,nz,4)
c      call wtap1 (agaut,bnz,nz,4)
c      call wtap1 (asdep,bnz,nz,4)
c      call wtap1 (agsub,bnz,nz,4)
c      call wtap1 (apern,bnz,nz,4)
c      call wtap1 (aqls,bnz,nz,4)
c      call wtap1 (altrans,bnz,nz,4)
c      call wtap1 (actrans,bnz,nz,4)
c      call wtap1 (total_mq,bnz,nz,4)
c      call wtap1 (agwet,bnz,nz,4)
c      call wtap1 (afgr,bnz,nz,4)
c      call wtap1 (avelw,bnz,nz,4)
c      call wtap1 (avesw,bnz,nz,4)
c      call wtap1 (amelt,bnz,nz,4)
c      call wtap1 (aother,bnz,nz,4)
c      call wtap1 (total_mt,bnz,nz,4)
c      call wtap1 (tavet,bnz,nz,4)
c      call wtap1 (qavet,bnz,nz,4)
c      call wtap1 (t_sfcq,bnnt,nnt,4)
c      call wtap1 (t_sfct,bnnt,nnt,4)
c      call wtap1 (t_ub,bnnt,nnt,4)
c      call wtap1 (t_microq,bnnt,nnt,4)
c      call wtap1 (t_largeq,bnnt,nnt,4)
c      call wtap1 (t_larget,bnnt,nnt,4)
c      call wtap1 (t_microt,bnnt,nnt,4)
c      call wtap1 (tsfc_1,bnx,nx,4)
c      call wtap1 (qsfc_1,bnx,nx,4)
c      call wtap1 (tlsw,bnz,nz,4)
c      call wtap1 (qlsw,bnz,nz,4)
c      call wtap1 (ttlsw,bnz,nz,4)
c      call wtap1 (tqlsw,bnz,nz,4)
c      call wtap1 (tcof,bnz,nz,4)
c      call wtap1 (textra,bnz,nz,4)
c      call wtap1 (qextra,bnz,nz,4)
c      call wtap1 (tlsw1,bnz,nz,4)
c      call wtap1 (qlsw1,bnz,nz,4)
c      call wtap1 (ttlsw1,bnz,nz,4)
c      call wtap1 (tqlsw1,bnz,nz,4)
c2223  continue
      endif

C     Write to restart file

      call wtap (dpt,3)
      call wtap (dqv,3)
      call wtap (qcl,3)
      call wtap (qrn,3)
      call wtap (qci,3)
      call wtap (qcs,3)
      call wtap (qcg,3)
      call wtap (u,3)
      call wtap (v,3)
      call wtap (w,3)
      call wtap (ak,3)
      call wtap (dpt1,3)
      call wtap (dqv1,3)
      call wtap (qcl1,3)
      call wtap (qrn1,3)
      call wtap (qci1,3)
      call wtap (qcs1,3)
      call wtap (qcg1,3)
      call wtap (uu1,3)
      call wtap (vv1,3)
      call wtap (ww1,3)
      call wtap (ak1,3)
      call wtap (umd,3)
      call wtap (vmd,3)
      call wtap (wmd,3)
      call wtap (rsw,3)
      call wtap (rlw,3)
      call wtap (ppress,3)
      call wtapt (3)
        rdt=dt
       call wstart (3)
c      call wstrad (3,cosz,month,iday,hrl)
       call wstrad (3,cosz,month,iday,hrl,rlat)  ! 8/21/01 shie
       call wstrad2 (3)  ! newly added! 8/21/01 shie
       write(6,12345) n,isec,nran,sec,aminut,rdt
cli   call outrf ( ircl,npp1,pl,pa)
c     *************************************
       anv1=0.
       cnv1=0.
       anv2=0.
       cnv2=0.
       anv3=0.
       cnv3=0.
       anv4=0.
       cnv4=0.
      do 110 ix=1,nxi
      do 110 jt=1,nt
       a1=rst(jt,ix)
       rsv(jt,ix)=0.0
       if(a1.ge.1.e-3) rsv(jt,ix)=23.62+12.5*log10(a1)
       a2=rsv(jt,ix)
       if(a2.gt.29.) cnv1=cnv1+rst(jt,ix)
       if(a2.le.29.) anv1=anv1+rst(jt,ix)
       if(a2.gt.39.) cnv2=cnv2+rst(jt,ix)
       if(a2.le.39.) anv2=anv2+rst(jt,ix)
       if(a1.gt.10.) cnv3=cnv3+rst(jt,ix)
       if(a1.le.10.) anv3=anv3+rst(jt,ix)
       if(a1.gt.20.) cnv4=cnv4+rst(jt,ix)
       if(a1.le.20.) anv4=anv4+rst(jt,ix)
  110 continue
       write(6,1234) cnv1,anv1,cnv2,anv2,cnv3,anv3,cnv4,anv4
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,'  time = ',isec,'  seconds'
cccccccccccccccccccccccccccccccccccccccc
c      call dom_ave (aminut,tave,denq,thee)
cccccccccccccccccccccccccccccccccccccccc
c        d22t=d2t
c        if(ijkadv .eq. 1) d22t=dt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        dt_stat=ndt_stat
c        PRINT *,'                   '
c        PRINT *,'     WATER VAPOR BUDGET '
c        PRINT *,'                   '
c
c      call budget_m (dt_stat)
c
c      do k=2,kles
c        altrans(k)=am(k)*(sqf(k)-sqf(k+1))*rdz*rrho(k)
c        atsqf1(k)=tsqf1(k,1,4)
c        atsqf2(k)=tsqf2(k,1,4)
c        atsqf (k)=tsqf (k,1,4)
c        avefluxq(k)=am(k)*(tsqf(k,1,4)-tsqf(k+1,1,4))*rd2z*rrho(k)
c      end do

c        PRINT *,'TSQF1 (LARGE SCALE FORCING - advect)'
c      call var_ave(atsqf1,atsqf1_t)
c        PRINT *,'TSQF1 (LARGE SCALE FORCING) = ',atsqf1_t
c
c        PRINT *,'tsqf2 (LOCAL TIME CHANGE - advect)'
c      call var_ave(atsqf2,atsqf2_t)
c        PRINT *,'tsqf2 (LOCAL TIME CHANGE)   = ',atsqf2_t
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        PRINT *,'CONDENSATION'
c      call var_ave (acoc,t)
c        PRINT *,'EVAPORATION'
c      call var_ave (acoe,t)
c        PRINT *,'DEPOSITION OF QI'
c      call var_ave (agaut,t)
c        PRINT *,'DEPOSITION OF QS'
c      call var_ave (asdep,t)
c        PRINT *,'DEPOSITION OF QG'
c      call var_ave (agsub,t)
c        PRINT *,'SUBLIMATION OF QI'
c      call var_ave (apern,t)
c        PRINT *,'SUBLIMATION OF QG TO QV'
c      call var_ave (agwet,t)
c        PRINT *,'SUBLIMATION OF QS TO QV'
c      call var_ave (afgr,t)
c        PRINT *,'NET CONDENSATION'
c      call var_ave (total_mq,t)
c
c        PRINT *,'NET MELTING FROM QS AND QG'
c      call var_ave (amelt,t)
c        PRINT *,'NET FREEZING'
c      call var_ave (aother,t)
C
c      PRINT *,'WB DQV/DZ - advect (extra)  '
c      call var_ave (qextra,t)
C
c        PRINT *,'TOTAL LARGE-SCALE FORCINGS FROM MAIN PROGRAM'
c        PRINT *,' w d(qb)/dz '
c      call var_ave (tqlsw,t)
c
c        PRINT *,'wbar d(qb)/dz from main program'
c      call var_ave (tqlsw1,t)
c
c        PRINT *,'total forcings from main program'
c      call var_ave (qls,t)
C
c       PRINT *,'WB D(DQV)/DZ FROM SQF = '
c      call var_ave (altrans,t)
C
c       PRINT *,'WB D(DQV)/DZ FROM TSQF = '
c      call var_ave (avefluxq,t)
C
c        SQFTT=10.*SQF(2)
c      PRINT *,'LATENT HEAT FLUXES = ',SQFTT
C
C
cccccccccccccccccccccccccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        PRINT *,'                   '
c        PRINT *,'      HEAT BUDGET  '
c        PRINT *,'                   '
c
c      do k=2,kles
c        altrans(k)=am(k)*(stf(k)-stf(k+1))*rdz*rrho(k)
c        atsqf1(k)=tstf1(k,1,4)
c        atsqf2(k)=tstf2(k,1,4)
c        atsqf (k)=tstf (k,1,4)
c        avefluxq(k)=am(k)*(tstf(k,1,4)-tstf(k+1,1,4))*rd2z*rrho(k)
c      end do
c
c        stf(2)=stf(2)/denave


c        PRINT *,'TSTF1 (LARGE SCALE FORCING - advect)'
c      call var_avet(atsqf1,atstf1_t)
c        PRINT *,'TSTF1 (LARGE SCALE FORCING) = ',atstf1_t
c
c        PRINT *,'TSTF2 (LOCAL TIME CHANGE - advect)'
c      call var_avet(atsqf2,atstf2_t)
c        PRINT *,'TSTF2 (LOCAL TIME CHANGE)   = ',atstf2_t
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        PRINT *,'CONDENSATION'
c      call var_avet(acoc,t)
c        PRINT *,'EVAPORATION'
c      call var_avet(acoe,t)
c        PRINT *,'DEPOSITION OF QI'
c      call var_avet(agaut,t)
c        PRINT *,'DEPOSITION OF QS'
c      call var_avet(asdep,t)
c        PRINT *,'DEPOSITION OF QG'
c      call var_avet(agsub,t)
c        PRINT *,'SUBLIMATION OF QI'
c      call var_avet(apern,t)
c        PRINT *,'SUBLIMATION OF QG TO QV'
c      call var_avet(agwet,t)
c        PRINT *,'SUBLIMATION OF QS TO QV'
c      call var_avet(afgr,t)
c        PRINT *,'NET MELTING FROM QS AND QG'
c      call var_avet(amelt,t)
c        PRINT *,'NET FREEZING'
c      call var_avet(aother,t)
c
c        PRINT *,'NET LONGWAVE COOLING'
c      call var_avet(avelw,t)
c        PRINT *,'NET SHORTWAVE HEATING'
c      call var_avet(avesw,t)
c        PRINT *,'NET HEATING FROM MICROPHYSICAL PROCESSES'
c      call var_avet(total_mt,t)
c
c        PRINT *,'WB DPT/DZ 2DT (extra) - advect '
c      call var_avet(textra,t)
c
c        PRINT *,'TOTAL LARGE-SCALE FORCINGS FROM MAIN PROGRAM'
c        PRINT *,' w d(tb)/dz '
c      call var_avet(ttlsw,t)
c
c        PRINT *,'wbar d(tb)/dz from main program'
c      call var_avet(ttlsw1,t)
c
c        PRINT *,'total forcings from main program'
c      call var_avet(tls,t)

c        PRINT *,'TOTAL CORIOLIS FORCINGS FROM MAIN PROGRAM'
c      call var_avet(tcof,t)
C
c        PRINT *,'WB D(DPT)/DZ FROM STF = '
c      call var_avet(altrans,altrans_t)
C
c        PRINT *,'WB D(DPT)/DZ FROM TSQF = '
c      call var_avet(avefluxq,avefluxt_t)

c        StFTT=StF(2)
c      PRINT *,'SENSIBLE HEAT FLUXES = ',STFTT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      call budget (taa)
      write(6,1099)
      stop
  201 format(10h time step,i5,6h time=,f8.2,9h min smt=,e11.4,5h smv=,
     1  e11.4,5h smq=,e11.4,5h smc=,e11.4,5h sme=,e11.4,5h smr=,e11.4)
  202 format(10h time step,i5,6h time=,f8.2,9h min sqc=,e11.4,5h sqr=,
     1  e11.4,5h sqi=,e11.4,5h sqs=,e11.4,5h sqg=,e11.4)
  207 FORMAT(10H TIME STEP,I5,6H TIME=,F8.2,10H MIN SMRC=,E11.4)
  125 format(1x,2hlv,6x,3h st,6x,3h sv,6x,3h sq,6x,3h sc,6x,3h se,6x,
     1   3hsqa,6x,3hcld,5x,4hcldu,5x,4hcldd,5x,3h ub,6x,3h vb,6x,3h wb)
 1025 format(1x,2hlv,6x,3hQ1Z,6x,3hQ2Z,6x,3hQ1H,6x,3hQ2H)
  126 format(1x,2hlv,6x,3hsqc,6x,3hsqr,6x,3hsqi,6x,3hsqs,6x,3hsqg,6x,
     1   3htqc,6x,3htqr,6x,3htqi,6x,3htqs,6x,3htqg)
  133 format(1x,4h l= ,i4,3x,2hj=,i3,6x,4hcon=,i4,6x,4hanv=,i4,6x,i8)
  134 format(1x,i2,1p14e9.2)
  135 format(129h lv     ihom      idw     imlt     saut     saci     sa
     1cw     raci     iacr     raut     racw      sfw      sfi      gacs
     2    gacw)
  136 format(111h lv     gaci     gacr     gwet     gaut     racs     sa
     1cr      gfr     smlt     gmlt     sdep     ssub     gsub)
  137 format(57h lv     subi     d3ri     d3ir     d2sr     depi     gdr
     1y,5x,4herns,5x,4hwgrs,5x,4hqsws)
  138 format(129h lv      coc      coe      smf      qc0      qr0      q
     1i0      qs0      qg0     sqc0     sqr0     sqi0     sqs0      sqg0
     2        )
  139 format(129h lv      smf      smu      smd      stf      stu      s
     1td      sqf      squ      sqd      sdt      sdq      stv      sft
     2     sfq            )
  140 format(129h lv      stf     stf2     stf1      sqf     sqf2     sq
     1f1      ut1      ut2      vt1      vt2                            
     2                    )
  154 format(6h nrun=,i11/
     1  6h il2 =,i11/6h kl2 =,i11/4h dx=,-2pf13.2/4h dz=,f13.2/
     2  6h time=,0pf11.2/4h dt=,f13.2/6h psfc=,-3pf11.2/5h rd1=,
     3  3pf12.2/5h rd2=,0pf12.2/7h bound=,3pf10.2/4h ck=,0pf13.2/
     4  4h ce=,f13.2/5h rm0=,f12.2/5h rm1=,f12.2/5h eps=,f12.2/
     5  5h inp=,i12/5h imp=,i12/6h iwp1=,i11/
     6  6h ised=,i11/6h itrn=,i11/6h ngtp=,i11/6h coto=,3pf11.2)
  155 format(//42h two-d slab-symmetric cloud ensemble model///
     1  23h second order flux form///)
  160 format(/16h ri = mm per hr ,/)
  203 format(1x,i2,128a1)
  204 format(1x,i4,40i3)
 1099 format(////5x,15hend of this run)
 1234 format(////5x,1p8e10.3)
 2345 format(/3x,8i10)
 3456 format(3x,128i1)
12345 format(////5x,3i9,3x,3f12.2)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advect (x,x1,aa,AA1,ism,ist,isq,iss,isr,isuv,dlt1,dlt2)
      parameter (NX=514,NZ=43,ITT=244,NT=2880)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON/SFLUXS/ SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      common/itoga/ itoga,ISFC,ICE,ICE2
      common/cpsbm/ icps_bm,iexplicit
      common/timestat/ ndt_stat,itime_ave,mmave
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      common/bsfc/ tsfc_1(nx),qsfc_1(nx)
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bcor/ irf,iadvh,irfg,idq,ismg
      common/bcorr/ id
      common/bstart/ n,isec,nran
      common/bsize/ it1,it2,kt1,kt2

      common/add_f/ i_four

      common/o4/ a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
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
      common/b1u/ u(nx,nz)
      common/b1v/ v(nx,nz)
      common/b1w/ w(nx,nz)
      common/b1a/ ak(nx,nz)
      common/b2u/ uu1(nx,nz)
      common/b2v/ vv1(nx,nz)
      common/b2w/ ww1(nx,nz)
      common/bsat/ xxw(nx,nz)
      common/bsat1/ aaa(nx,nz)
      common/badv/ dfc(nx,nz)
      common/btv/ vtp(nx,nz)
      common/dumuw/ umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
cccccc
c      common/bstsi/ ceds1i(nx,nz,4),ceds2i(nx,nz,4),
c     1   tstfi(nx,nz,4),tsqfi(nx,nz,4),rsli(nx,nz,4)
      common/bw/ xxx(nx,nz),dfx(nx,nz),trah(nx,nz),trav(nx,nz)
      common/bw00/ y(nx,nz),xx0(nx,nz),xx00(nx,nz),xx000(nx,nz)
cccccc
      common/b4/ tbsk(nz),bsk(nz),bsk2(nz),bskm(nz),bskm4(nz),
     1        bsi(nz),bsi2(nz),bsim(nz),bsim4(nz)
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/b8/ smf(nz),smu(nz),smd(nz),smn(nz),stf(nz),stu(nz),
     $   std(nz),stn(nz),sqf(nz),squ(nz),sqd(nz),sqn(nz),sdt(nz),
     $   sdq(nz),stv(nz),
     $   stf1(nz),stu1(nz),std1(nz),stn1(nz),sqf1(nz),squ1(nz),
     $   sqd1(nz),sqn1(nz),stf2(nz),stu2(nz),std2(nz),stn2(nz),sqf2(nz),
     $   squ2(nz),sqd2(nz),sqn2(nz),cld(nz),cldu(nz),cldd(nz),cldn(nz)
      common/b9/ sut1(nz),suc1(nz),sun1(nz),suu1(nz),sud1(nz),aub(nz),
     1   ssu(nz),sut2(nz),suc2(nz),sun2(nz),suu2(nz),sud2(nz),auc1(nz),
     2   aun1(nz),auu1(nz),aud1(nz),auc2(nz),aun2(nz),auu2(nz),aud2(nz),
     3   svt1(nz),svc1(nz),svn1(nz),svu1(nz),svd1(nz),avb(nz),
     4   ssv(nz),svt2(nz),svc2(nz),svn2(nz),svu2(nz),svd2(nz),avc1(nz),
     5   avn1(nz),avu1(nz),avd1(nz),avc2(nz),avn2(nz),avu2(nz),avd2(nz)

      COMMON/BTT/ S1(16,NZ),S2(16,NZ),S3(16,NZ),S4(16,NZ),S5(16,NZ),
     1   S6(16,NZ),S7(16,NZ),S8(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     2   S14(16,NZ),S15(16,NZ),S16(16,NZ),SN1(5,NZ),SN2(5,NZ),SN3(5,NZ),
     3   SN4(5,NZ),SN5(5,NZ),SN6(5,NZ),SN7(5,NZ),SN8(5,NZ),SN11(5,NZ),
     4   SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),SN15(5,NZ),SN16(5,NZ)


      common/bsts1/ tut1(nz,4,7),tut2(nz,4,7),tvt1(nz,4,7),tvt2(nz,4,7),
     $ tstf(nz,4,7),tstf1(nz,4,7),tstf2(nz,4,7),tsqf(nz,4,7),
     $ tsqf1(nz,4,7),tsqf2(nz,4,7),tsqq(nz,4,7),tsqq1(nz,4,7),
     $ fcld(nz,4,7)
      common/bsts2/ snth(nz,4,7),sntv(nz,4,7),sntd(nz,4,7),snqh(nz,4,7),
     $        snqv(nz,4,7),snqd(nz,4,7),snhh(nz,4,7),snhv(nz,4,7),
     $        snhd(nz,4,7)
      common/bsts20/ othert_add(nz,4,7),otherq_add(nz,4,7)
      common/qcbf/ qcl_buffer(nx,nz)
      common/gbs/ tlsw(nz),qlsw(nz),ttlsw(nz),tqlsw(nz),smfff(nz),
     1   smffff(nz),tcof(nz),textra(nz),qextra(nz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/ba/y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),y8(nx)
     $   ,y9(nx),y10(nx),y11(nx),y12(nx),y13(nx),y14(nx),y15(nx),y16(nx) 
      common/bls/ y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/damp/ rfa(nz),rfa1(nz),cntc(nz),cgwd
      common/bch/ it(nx),iv(nt),ics(nx,4),ibz(nx,4)
      common/surface/ sun_4(nx,4)
      common/bch1/ rby(7)
c-----------------------------------------------------------------------
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
      real    x(nx,nz),x1(nx,nz),aa(nz),AA1(NZ),y17(nx),y18(nx)
      real    y19(nx),y20(nx),y21(nx),y22(nx),y23(nx)
      save
c     ***   statistics for cloud properties   ******
      BARFCT=2.*100.*RHO1(2)    !convert from MKS to cgs
       d22t=d2t
       assss=1.
       awks=1.
       d2t=d2t
       if (ijkadv .eq. 1) then
         d2t=dt
         assss=0.
         awks=.5
       endif

       a0000=0.
         if (IMLIFTING .eq. 1) a0000=1.
       a_irf=0.
         if (irf .eq. 1) a_irf=1.
c
      do k=1,kmax
      do i=1,imax
        trah(i,k)=0.0
        trav(i,k)=0.0
        y(i,k)=0.0
        xx0(i,k)=0.0
        xx00(i,k)=0.0
        xx000(i,k)=0.0
      enddo
      enddo
c
      do 11 i=1,imax
       y19(i)=0.0
       y20(i)=0.0
       y21(i)=0.0
       y22(i)=0.0
       y23(i)=0.0
       y17(i)=0.0
   11  y18(i)=0.0

      A0 =ndt_stat*0.5*RIL2
      A20=ndt_stat*    RIL2
      a300=d2t*ril2
      if(ism+ist+isq+idq.eq.0) go to 3333
      y1(2)=0.
      y1(kmax)=0.
      a0=ndt_stat*0.5*ril2
      if (idq.eq.0) then
      do 12 k=2,kles
       y15(k)=ndt_stat*ril2*rdz*am1(k)*rrho(k)
       y7(k)=.5*(ub(k)+ub(k-1))
   12  y8(k)=.5*(vb(k)+vb(k-1))
      do 600 i=2,iles
      do 60 k=2,kles
       y1(k)=rho1(k)*(ww1(i,k)+wb(k)*y0(i))
       y2(k)=rho1(k)*ww1(i,k)
       y3(k)=.5*(x(i,k)+x(i,k-1))
       y4(k)=y1(k)*y3(k)
       y5(k)=.25*(uu1(i,k)+uu1(i+1,k)+uu1(i,k-1)+uu1(i+1,k-1))
       y6(k)=.25*(vv1(i,k)+vv1(i+1,k)+vv1(i,k-1)+vv1(i+1,k-1))
       y9(k)=y5(k)-y7(k)
       y10(k)=y6(k)-y8(k)
       y11(k)=y1(k)*y9(k)
       y12(k)=y1(k)*y10(k)
       y13(k)=y1(k)*y5(k)
       y14(k)=y1(k)*y6(k)
   60 continue
      do 61 k=2,kles
        y17(k)=qcl(i,k)+qci(i,k)
        y18(k)=qcl(i,k-1)+qci(i,k-1)
       if(y17(k).lt.1.e-5 .or. y18(k).lt.1.e-5) then
         y16(k)=rho1(k)*(2.*assss*(ww1(i,k)+wb(k)*y0(i))*y3(k)
     1          -((ak(i,k)*(1.+xxw(i,k))+ak(i,k-1)*(1.+xxw(i,k-1)))
     2          *am1(k)*(x1(i,k)-x1(i,k-1)+aa(k)-aa(k-1))
     3          +awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,k-1)))*rdz)*a0
       else
         y16(k)=rho1(k)*(2.*assss*(ww1(i,k)+wb(k)*y0(i))*y3(k)+dfc(i,k)
     1          -awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,k-1))*rdz)*a0
       endif
   61 continue
      do 30 k=2,kles
      qq=qcl(i,k)+qrn(i,k)+qcl(i,k-1)+qrn(i,k-1)
     1   +qci(i,k)+qcs(i,k)+qcg(i,k)+qci(i,k-1)+qcs(i,k-1)+qcg(i,k-1)
ccccc
         WW=2.*A0*Y1(K)
        iww=0
        jww=0
        rq=.01*(ww1(i,k)+wb(k))
         if (rq .gt. -0.5) then
          iww=min(rq+0.5, 15.) + 1
         else
          jww=max(rq-0.5, -5.0)
          jww=abs(jww)
         endif
ccccc
      if(ism.eq.0) go to 10
       smf(k)=smf(k)+ww
       sut1(k)=sut1(k)+y11(k)
       sut2(k)=sut2(k)+y13(k)
       svt1(k)=svt1(k)+y12(k)
       svt2(k)=svt2(k)+y14(k)
      if(qq.lt.2.e-5) go to 601
       cld(k)=cld(k)+1.
       suc1(k)=suc1(k)+y11(k)
       suc2(k)=suc2(k)+y13(k)
       auc1(k)=auc1(k)+y5(k)
       auc2(k)=auc2(k)+y9(k)
       svc1(k)=svc1(k)+y12(k)
       svc2(k)=svc2(k)+y14(k)
       avc1(k)=avc1(k)+y6(k)
       avc2(k)=avc2(k)+y10(k)
      if(ww.lt.0.) go to 602
       cldu(k)=cldu(k)+1.
       smu(k)=smu(k)+ww
       suu1(k)=suu1(k)+y11(k)
       suu2(k)=suu2(k)+y13(k)
       auu1(k)=auu1(k)+y5(k)
       auu2(k)=auu2(k)+y9(k)
       svu1(k)=svu1(k)+y12(k)
       svu2(k)=svu2(k)+y14(k)
       avu1(k)=avu1(k)+y6(k)
       avu2(k)=avu2(k)+y10(k)
      if(iww.eq.0) go to 603
       s1(iww,k)=s1(iww,k)+ww
       s2(iww,k)=s2(iww,k)+1.
       s11(iww,k)=s11(iww,k)+y9(k)
       s12(iww,k)=s12(iww,k)+y10(k)
       s13(iww,k)=s13(iww,k)+y5(k)
       s14(iww,k)=s14(iww,k)+y6(k)
       s15(iww,k)=s15(iww,k)+y11(k)
       s16(iww,k)=s16(iww,k)+y12(k)
      go to 603
  602  cldd(k)=cldd(k)+1.
       smd(k)=smd(k)+ww
       sud1(k)=sud1(k)+y11(k)
       sud2(k)=sud2(k)+y13(k)
       aud1(k)=aud1(k)+y5(k)
       aud2(k)=aud2(k)+y9(k)
       svd1(k)=svd1(k)+y12(k)
       svd2(k)=svd2(k)+y14(k)
       avd1(k)=avd1(k)+y6(k)
       avd2(k)=avd2(k)+y10(k)
      if(jww.eq.0) go to 603
       sn1(jww,k)=sn1(jww,k)+ww
       sn2(jww,k)=sn2(jww,k)+1.
       sn11(jww,k)=sn11(jww,k)+y9(k)
       sn12(jww,k)=sn12(jww,k)+y10(k)
       sn13(jww,k)=sn13(jww,k)+y5(k)
       sn14(jww,k)=sn14(jww,k)+y6(k)
       sn15(jww,k)=sn15(jww,k)+y11(k)
       sn16(jww,k)=sn16(jww,k)+y12(k)
      go to 603
  601  cldn(k)=cldn(k)+1.
       smn(k)=smn(k)+ww
       sun1(k)=sun1(k)+y11(k)
       sun2(k)=sun2(k)+y13(k)
       aun1(k)=aun1(k)+y5(k)
       aun2(k)=aun2(k)+y9(k)
       svn1(k)=svn1(k)+y12(k)
       svn2(k)=svn2(k)+y14(k)
       avn1(k)=avn1(k)+y6(k)
       avn2(k)=avn2(k)+y10(k)
  603 continue
   10 if(ist.eq.0) go to 20
       stf1(k)=stf1(k)+y4(k)
       stf2(k)=stf2(k)+y3(k)
       stf(k)=stf(k)+y16(k)
       sdt(k)=sdt(k)+y1(k)*(aa(k)-aa(k-1))*y15(k)
      if(qq.lt.2.e-5) go to 611
      if(ww.lt.0.) go to 612
       stu1(k)=stu1(k)+y4(k)
       stu2(k)=stu2(k)+y3(k)
       stu(k)=stu(k)+y16(k)
      if(iww.eq.0) go to 613
       s3(iww,k)=s3(iww,k)+y4(k)
       s4(iww,k)=s4(iww,k)+y3(k)
       s5(iww,k)=s5(iww,k)+y16(k)
      go to 613
  612  std1(k)=std1(k)+y4(k)
       std2(k)=std2(k)+y3(k)
       std(k)=std(k)+y16(k)
       if(jww.eq.0) go to 613
       sn3(jww,k)=sn3(jww,k)+y4(k)
       sn4(jww,k)=sn4(jww,k)+y3(k)
       sn5(jww,k)=sn5(jww,k)+y16(k)
      go to 613
  611  stn1(k)=stn1(k)+y4(k)
       stn2(k)=stn2(k)+y3(k)
       stn(k)=stn(k)+y16(k)
  613 continue
      go to 30
   20 if(isq.eq.0) go to 30
       sqf1(k)=sqf1(k)+y4(k)
       sqf2(k)=sqf2(k)+y3(k)
       sqf(k)=sqf(k)+y16(k)
       sdq(k)=sdq(k)+y1(k)*(aa(k)-aa(k-1))*y15(k)
      if(qq.lt.2.e-5) go to 621
      if(ww.lt.0.) go to 622
       squ1(k)=squ1(k)+y4(k)
       squ2(k)=squ2(k)+y3(k)
       squ(k)=squ(k)+y16(k)
      if(iww.eq.0) go to 623
       s6(iww,k)=s6(iww,k)+y4(k)
       s7(iww,k)=s7(iww,k)+y3(k)
       s8(iww,k)=s8(iww,k)+y16(k)
      go to 623
  622  sqd1(k)=sqd1(k)+y4(k)
       sqd2(k)=sqd2(k)+y3(k)
       sqd(k)=sqd(k)+y16(k)
      if(jww.eq.0) go to 623
       sn6(jww,k)=sn6(jww,k)+y4(k)
       sn7(jww,k)=sn7(jww,k)+y3(k)
       sn8(jww,k)=sn8(jww,k)+y16(k)
      go to 623
  621  sqn1(k)=sqn1(k)+y4(k)
       sqn2(k)=sqn2(k)+y3(k)
       sqn(k)=sqn(k)+y16(k)
  623  continue
   30 continue
  600 continue
       do 9000 k=2,kles
        do 31 i=2,iles
   31   y20(i)=0.0
       if (k .eq. 2) then
        do 32 i=2,iles
cc          vvss =.5*abs( u(i,2)+u(i+1,2) )
cc          y4(i)=max(vvss, 400.E0)
c         y4(i)=      vvss
cc          y1(i)=1.e-3*rho1(2)*(1.1+.04*(.01*y4(i)))*y4(i)
cc        if (ist.eq.1)  y20(i)=max(y1(i)*(ts(i)-aa(2)-x1(i,2)) ,0.)
cc        if (isq.eq.1)  y20(i)=max(y1(i)*(qs(i)-aa(2)-x1(i,2)) ,0.)
c
        IF (ISmg.EQ.1)  Y20(I)=tsfc_1(i)
        IF (ISmg.EQ.2)  Y20(I)=qsfc_1(i)
        IF (ITOGA .EQ. 1) THEN
          if (ismg .eq. 1) y20(i)=100.*rho1(2)*swt(i)/pi(2)
          if (ismg .eq. 2) y20(i)=100.*rho1(2)*swq(i)
        ENDIF
   32   continue
       end if
        do 91 i=2,iles
        y1(i)=rho1(k)*(ww1(i,k)+wb(k)*y0(i))
        y2(i)=rho1(k)*ww1(i,k)
        y3(i)=.5*(x(i,k)+x(i,k-1))
        y4(i)=y1(i)*y3(i)
        y5(i)=.25*(uu1(i,k)+uu1(i+1,k)+uu1(i,k-1)+uu1(i+1,k-1))
        y6(i)=.25*(vv1(i,k)+vv1(i+1,k)+vv1(i,k-1)+vv1(i+1,k-1))
        y9(i)=y5(i)-y7(k)
        y10(i)=y6(i)-y8(k)
        y11(i)=y1(i)*y9(i)
        y12(i)=y1(i)*y10(i)
        y13(i)=y1(i)*y5(i)
        y14(i)=y1(i)*y6(i)
         y17(i)=qcl(i,k)+qci(i,k)
         y18(i)=qcl(i,k-1)+qci(i,k-1)
         y19(i)=0.
        y23(i)=y17(i)+y18(i)+qrn(i,k  )+qcs(i,  k)+qcg(i,  k)+
     1                       qrn(i,k-1)+qcs(i,k-1)+qcg(i,k-1)
ctao
        IF(IMLIFTING .EQ. 1) Y19(I)=-a20*AM(K)*Y0(I)*
     1                  (WB(K+1)+WB(K))*(X(I,K+1)-X(I,K-1))*RD4Z
        IF (IRF.EQ.1) Y19(I)=Y19(I)-a20*RFA(K)*X1(I,K)

cccccc  tao (11-13-97)   cccccccccccccccccccccccccccccccccccccccccccccc
        if (ijkadv .eq. 0) then
c          if(k.gt.3 .and. k.lt.kl2) then
c            y19(i)=y19(i)-0.002*(x1(i,k+2)+x1(i,k-2)-4.*(x1(i,k+1)
c     1                  +x1(i,k-1))+6.*x1(i,k))/dt
c          endif
c          if(k.eq.3) y19(i)=y19(i)+.008*(x1(i,4)-2.*x1(i,3)
c     1                           +x1(i,2))/dt
c          if(k.eq.kl2) y19(i)=y19(i)+.008*(x1(i,kles)-2.*x1(i,kl2)
c     1                      +x1(i,kl2-1))/dt
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctao
        if(y17(i).lt.1.e-5 .or. y18(i).lt.1.e-5) then
c         y16(i)=rho1(k)* (2.*assss*(ww1(i,k)+wb(k)*y0(i))*y3(i)
          y16(i)=rho1(k)* (2.*assss*ww1(i,k)*y3(i)
     1          -((ak(i,k)*(1.+xxw(i,k))+ak(i,k-1)*(1.+xxw(i,k-1)))
     2          *am1(k)*(x1(i,k)-x1(i,k-1)+aa(k)-aa(k-1))
     3          +awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,k-1)))*rdz)*a0
     4          +y20(i)*a20
          
        else
c         y16(i)=rho1(k)*(2.*assss*(ww1(i,k)+wb(k)*y0(i))*y3(i)+dfc(i,k)
          y16(i)=rho1(k)*(2.*assss*ww1(i,k)*y3(i)+dfc(i,k)
     1          -awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,k-1))*rdz)*a0
     2          +y20(i)*a20
        endif
          y15(i)=.01*(ww1(i,k)+wb(k)*y0(i))
        if (ist.eq.1) then
          y21(i)= (x(i,k)+ta(k)-x1(i,k)-ta1(k))/dt
C         y22(i)=(am1(k)*fd(k)*(ww1(i,k)+wb(k)*y0(i))
C    1     +am1(k+1)*fd(k+1)*(ww1(i,k+1)+wb(k+1)*y0(i)))*rdz*rrho(k)
          y22(i)=(am1(k)*fd(k)*ww1(i,k)
     1           +am1(k+1)*fd(k+1)*ww1(i,k+1))*rdz*rrho(k)
        end if
        if (isq.eq.1) then
          y21(i)= (x(i,k)+qa(k)-x1(i,k)-qa1(k))/dt
C          y22(i)=(am1(k)*fe(k)*(ww1(i,k)+wb(k)*y0(i))
C    1     +am1(k+1)*fe(k+1)*(ww1(i,k+1)+wb(k+1)*y0(i)))*rdz*rrho(k)
          y22(i)=(am1(k)*fe(k)*ww1(i,k)
     1           +am1(k+1)*fe(k+1)*ww1(i,k+1))*rdz*rrho(k)
        end if
   91  continue




       do 900 kc=1,7
         do 95 mt=1,4
         do 95 i=2,iles
          ibz(i,mt)=0
          ibz(i,1)=1
           if(ics(i,mt).eq.1) ibz(i,mt)=1
           if(kc.eq.4) go to 93
           if(kc.le.3) go to 96
            if (y15(i).gt.rby(kc)) ibz(i,mt)=0
            go to 93
   96       if (y15(i).lt.rby(kc)) ibz(i,mt)=0
   93     continue
   95     continue
         do 98 mt=1,4
          a1=0.
          a2=0.
          a3=0.
          a4=0.
          a5=0.
          a6=0.
          a7=0.
          a8=0.
          a9=0.
          a10=0.
          a45=0.
          a55=0.
          a66=0.
         do 99 i=2,iles
           if(ism.eq.1 .and. ibz(i,mt).eq.1) then
             a1=a1+y11(i)
             a2=a2+y13(i)
             a3=a3+y12(i)
             a4=a4+y14(i)

            if(y23(i) .ge. 2.e-5) a45=a45+1.
           endif
           if(ist.eq.1 .and. ibz(i,mt).eq.1) then
             a5=a5+y16(i)
c            if(kc.eq.4)then
c             tstfi(i,k,mt)=tstfi(i,k,mt)+y16(i)    
c            endif
c            a6=a6+y3(i)
             a6=a6+y21(i)*a20
c            a7=a7+y4(i)
             a7=a7+y22(i)*a20
             A55=A55+Y19(I)
           endif
           if(isq.eq.1 .and. ibz(i,mt).eq.1) then
             a8=a8+y16(i)
c            if(kc.eq.4)then
c             tsqfi(i,k,mt)=tsqfi(i,k,mt)+y16(i)      
c            endif
c            a9=a9+y3(i)
             a9=a9+y21(i)*a20
c            a10=a10+y4(i)
             a10=a10+y22(i)*a20
             A66=A66+Y19(I)
           endif
   99    continue
           tut1(k,mt,kc)=tut1(k,mt,kc)+a1
           tut2(k,mt,kc)=tut2(k,mt,kc)+a2
           tvt1(k,mt,kc)=tvt1(k,mt,kc)+a3
           tvt2(k,mt,kc)=tvt2(k,mt,kc)+a4

           othert_add(k,mt,kc)=othert_add(k,mt,kc)+a55
           otherq_add(k,mt,kc)=otherq_add(k,mt,kc)+a66

           tstf(k,mt,kc)=tstf(k,mt,kc)+a5
           tstf2(k,mt,kc)=tstf2(k,mt,kc)+a6
           tstf1(k,mt,kc)=tstf1(k,mt,kc)+a7

           tsqf(k,mt,kc)=tsqf(k,mt,kc)+a8
           tsqf2(k,mt,kc)=tsqf2(k,mt,kc)+a9
           tsqf1(k,mt,kc)=tsqf1(k,mt,kc)+a10

            fcld(k,mt,kc)= fcld(k,mt,kc)+a45
   98    continue
  900   continue
 9000   continue
       endif
c     ****** statistic for u- and v-momentum
      if(isuv.eq.1) then
       do 70 k=2,kles
        aub(k)=aub(k)+y7(k)
        avb(k)=avb(k)+y8(k)
   70  continue
      endif
c     ****** statistic for w(qc+qr+qi+qs+qg)
      if (idq.eq.1) then
      do 7000 k=2,kles
       do 71 i=2,iles
        y1(i)=rho1(k)*ww1(i,k)*(x(i,k)+x(i,k-1))*a0
        y2(i)=.01*ww1(i,k)
        y3(i)=rho1(k)*(assss*ww1(i,k)*(x(i,k)+x(i,k-1))
     1                 -(am1(k)*(x1(i,k)-x1(i,k-1))
     2          *(ak(i,k)*(1.+xxw(i,k))+ak(i,k-1)*(1.+xxw(i,k-1)))
     3            +awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,k-1)))*rdz)*a0
        if(isr.ne.1) then

          if (ismg .ne. 4) then
            y17(i)=dlt1*qcl(i,k)+dlt2*qci(i,k)
            y18(i)=dlt1*qcl(i,k-1)+dlt2*qci(i,k-1)
          else
            y17(i)=dlt1*qcl_buffer(i,k)+dlt2*qci(i,k)
            y18(i)=dlt1*qcl_buffer(i,k-1)+dlt2*qci(i,k-1)
          end if

         if (y17(i) .ge. 1.e-5 .and. y18(i) .ge. 1.e-5) then
          y3(i)=rho1(k)*(assss*ww1(i,k)*(x(i,k)+x(i,k-1))+dfc(i,k)
     1          -awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,k-1))*rdz)*a0
         endif
        endif
          y21(i)= (x(i,k)-x1(i,k))/dt
   71  continue
        do 700 kc=1,7
         do 75 mt=1,4
         do 75 i=2,iles
          ibz(i,mt)=0
          ibz(i,1)=1
           if(ics(i,mt).eq.1) ibz(i,mt)=1
           if(kc.eq.4) go to 73
           if(kc.le.3) go to 76
            if (y2(i).gt.rby(kc)) ibz(i,mt)=0
            go to 73
   76       if (y2(i).lt.rby(kc)) ibz(i,mt)=0
   73     continue
   75     continue
         do 89 mt=1,4
          a11=0.
          a12=0.
         do 79 i=2,iles
          if(ibz(i,mt).eq.1) then
c            a11=a11+y1(i)
             a11=a11+y21(i)*a20
             a12=a12+y3(i)
          endif
   79   continue
        tsqq(k,mt,kc)=tsqq(k,mt,kc)+a11
        tsqq1(k,mt,kc)=tsqq1(k,mt,kc)+a12
   89   continue
  700  continue
 7000  continue
      endif
c     ******   vertical advection terms   **************************
 3333 continue
c
cc
      if (ismg .le. 2) then
        do k=2,kles
        do i=2,iles
            aextra=am(k)*y0(i)*(wb(k+1)+wb(k))*
     1                          (x(i,k+1)-x(i,k-1))*rd4z*a300
          if (ismg .eq. 1) textra(k)=textra(k)-aextra
          if (ismg .eq. 2) qextra(k)=qextra(k)-aextra
        enddo
        enddo
      endif
cc
C-----put surface layer heat or moisture fluxs in
C
      DO I=2,ILES
        Y1(I)=0.
        Y3(I)=0.
      ENDDO
c
      IF (ISFC .EQ. 1 .AND. ITOGA .EQ. 1) THEN
        do i=2,iles
          if (ismg .eq. 1) then
            y3(i)=swt(i)*barfct/pi(2)
            SUN_4(I,3)=100.*swt(i)/pi(2)*.001
          endif
          if (ismg .eq. 2) then
            y3(i)=swq(i)*barfct
            SUN_4(I,4)=100.*swq(i)*.001
          endif
        enddo
      ENDIF
cc
      if (ijkadv .eq. 0) then
cc
c
      do 101 i=2,iles
  101  y1(i)=rho1(2)*ww1(i,2)*(x(i,2)+x(i,1))
      do 1000 k=3,kmax
        km=k-1
       do i=2,iles
         y2(i)=rho1(k)*ww1(i,k)*(x(i,k)+x(i,km))
       enddo

c       if(iss.eq.1) then
c       do 120 i=2,iles
c        y4(i)=rho1(k)*(-(am1(k)*(x1(i,k)-x1(i,km))
c     1        *(ak(i,k)*(1.+xxw(i,k))+ak(i,km)*(1.+xxw(i,km)))
c     2        +am1(k)*bsk2(k)*(x1(i,k)-x1(i,km)))*rdz)
c  120  continue
c       else
       do 140 i=2,iles
       y4(i)=rho1(k)*(
     1          -(am1(k)*(x1(i,k)-x1(i,km)+aa(k)-aa(km))*(ak(i,k)
     2          *(1.+xxw(i,k))+ak(i,km)*(1.+xxw(i,km)))
     3          +am1(k)*bsk2(k)*(x1(i,k)-x1(i,km)))*rdz)
  140  continue
c       endif

       if(isr.eq.0) then
       do 130 i=2,iles
         if (ismg .ne. 4) then
           y17(i)=dlt1*qcl(i,k)+dlt2*qci(i,k)
           y18(i)=dlt1*qcl(i,km)+dlt2*qci(i,km)
         else
           y17(i)=dlt1*qcl_buffer(i,k)+dlt2*qci(i,k)
           y18(i)=dlt1*qcl_buffer(i,km)+dlt2*qci(i,km)
         endif
        if (y17(i) .lt. 1.e-5 .or. y18(i) .lt. 1.e-5) go to 130
        y4(i)=rho1(k)*(dfc(i,k)-am1(k)*bsk2(k)*(x1(i,k)-x1(i,km))*rdz)
  130  continue
       endif
       if (km .eq. kles) then 
          do i=2,iles
             y2(i)=0.
             y4(i)=0.
          enddo
       endif
       do 170 i=2,iles
        y(i,km)=am(km)*(y1(i)-y2(i)+y3(i)-y4(i))*rd2z*rrho(km)
c         if(km .eq. kles) y(i,km)=am(km)*(y1(i)+y3(i))*rd2z*rrho(km)
c
c
c
c          if (ismg .eq. 1) then
c            q1_d_v(i,km)=q1_d_v(i,km)+am(km)*(y3(i)-y4(i))*rrho(km)*rd2z
c            q1_g_v(i,km)=q1_g_v(i,km)+am(km)*(y1(i)-y2(i))*rrho(km)*rd2z
c          endif
c
c          if (ismg .eq. 2) then
c            q2_d_v(i,km)=q2_d_v(i,km)+am(km)*(y3(i)-y4(i))*rrho(km)*rd2z
c            q2_g_v(i,km)=q2_g_v(i,km)+am(km)*(y1(i)-y2(i))*rrho(km)*rd2z
c          endif
c
c
          y1(i)=y2(i)
          y3(i)=y4(i)
  170  continue
 1000 continue
c        a2=am(kles)*rd2z*rrho(kles)
c       do i=2,iles
c          y(i,kles)=a2*(y2(i)+y4(i))
c       enddo

c     ******   horizontal advection terms   **************************

      do 2000 k=2,kles

         if(iadvh .eq. 2) then

            do i=2,imax
               y1(i)=uu1(i,k)*(x(i,k)+x(i-1,k))*rd2x
            enddo

            DO I=2,ILES
c
c               if (ismg .eq. 1) then
c                  q1_g_h(i,k)=q1_g_h(i,k)+(y1(i)-y1(i+1))
c               endif
c
c               if (ismg .eq. 2) then
c                  q2_g_h(i,k)=q2_g_h(i,k)+(y1(i)-y1(i+1))
c               endif
c
               Y(I,K)=Y(I,K)+Y1(I)-Y1(I+1)
            ENDDO

         else

c     ******   4-th order  horizontal advection terms   ***************

            y1(1)=uu1(2,k)*(d58x*(x(1,k)+x(2,k))-d16x*(x(3,k)+x(il2,k)))
 
            do i=2,il2
               y1(i)=uu1(i+1,k)*(d58x*(x(i,k)+x(i+1,k))
     1              -d16x*(x(i-1,k)+x(i+2,k)))
            enddo

            y3(1)=-d48x*uu1(1,k)*(x(1,k)+x(il2,k))

            do i=2,imax
               y3(i)=-d48x*uu1(i,k)*(x(i,k)+x(i-1,k))
            enddo

            do i=2,il2
               y(i,k)=y(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)

c               IF (ISMG .EQ. 1) THEN
c                  Q1_G_H(I,K)=Q1_G_H(I,K)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)
c               ENDIF
c               IF (ISMG .EQ. 2) THEN
c                  Q2_G_H(I,K)=Q2_G_H(I,K)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)
c               ENDIF

            enddo

            y(iles,k)=y(iles,k)+y1(il2)-y1(1)+y3(il2)-y3(3)

c            IF (ISMG .EQ. 1) THEN
c               Q1_G_H(ILES,K)=Q1_G_H(ILES,K)+y1(il2)-y1(1)+y3(il2)-y3(3)
c            ENDIF
c            IF (ISMG .EQ. 2) THEN
c               Q2_G_H(ILES,K)=Q2_G_H(ILES,K)+y1(il2)-y1(1)+y3(il2)-y3(3)
c            ENDIF

       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            do i=2,imax
               y1(i)=-(ak(i-1,k)+ak(i,k))*(x1(i,k)-x1(i-1,k))*r2dx2
            enddo

            do i=2,iles
               y(i,k)=y(i,k)+y1(i)-y1(i+1)
c
c
c            if (ismg .eq. 1) then
c               q1_d_h(i,k)=q1_d_h(i,k)+(y1(i)-y1(i+1))
c            endif
c
c            if (ismg .eq. 2) then
c               q2_d_h(i,k)=q2_d_h(i,k)+(y1(i)-y1(i+1))
c            endif
c
c
            enddo

            y1(1)=a4h(k)*(x1(3,k)-3.*(x1(2,k)-x1(1,k))-x1(il2,k))

            do i=2,il2
               y1(i)=a4h(k)*(x1(i+2,k)-3.*(x1(i+1,k)-x1(i,k))-x1(i-1,k))
            enddo

            do i=2,il2
c
C
c            IF (ISMG .EQ. 1) THEN
c               Q1_D_H(I,K)=Q1_D_H(I,K)+(Y1(i-1)-Y1(I))
c            ENDIF
c 
c            IF (ISMG .EQ. 2) THEN
c               Q2_D_H(I,K)=Q2_D_H(I,K)+(Y1(I-1)-Y1(I))
c            ENDIF
c
c
               y(i,k)=y(i,k)+y1(i-1)-y1(i)
            enddo

            y(iles,k)=y(iles,k)+y1(il2)-y1(1)

c            IF (ISMG .EQ. 1) THEN
c               Q1_D_H(ILES,K)=Q1_D_H(ILES,K)+(Y1(il2)-Y1(1))
c            ENDIF
c 
c            IF (ISMG .EQ. 2) THEN
c               Q2_D_H(ILES,K)=Q2_D_H(ILES,K)+(Y1(il2)-Y1(1))
c            ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DO 275 I=2,ILES
           Y(I,K)=Y(I,K)-a0000*(AM(K)*y0(i)*
     1                    (WB(K+1)+WB(K))*(X(I,K+1)-X(I,K-1))*RD4Z)
     2                  -a_irf*rfa(k)*x1(i,k)
  275   CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 2000 continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-13-97)
        IF (I_FOUR .EQ. 1) THEN
c
               r2dt2=0.005/dt
c
          do i=2,iles
             do k=2,kles
                y1(k)=x1(i,k+1)-2.*x1(i,k)+x1(i,k-1)
             enddo
             y1(1)=y1(2)
             y1(kmax)=y1(kles)
             do k=2,kles
                y(i,k)=y(i,k)-r2dt2*(y1(k+1)-2.*y1(k)+y1(k-1))
             enddo
          enddo
        ENDIF
C
c              IF (ISMG .EQ. 1) THEN
c                 do k=4,kl2-1
c                    Q1_D_V(I,K)=Q1_D_V(I,K)+r2dt2*(Y1(K-1)-Y1(K))
c                 enddo
c              ENDIF
c
c              IF (ISMG .EQ. 2) THEN
c                 do k=4,kl2-1
c                    Q2_D_V(I,K)=Q2_D_V(I,K)+r2dt2*(Y1(K-1)-Y1(K))
c                 enddo
c              ENDIF
c        IF (I_FOUR .EQ. 1) THEN
c
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c
c           do i=2,iles
c             do k=3,kl2-1
c               y1(k)=(x1(i,k+2)-3.*(x1(i,k+1)-x1(i,k))-x1(i,k-1))
c             enddo
c             do k=4,kl2-1
c                y(i,k)=y(i,k)+r2dt2*(y1(k-1)-y1(k))/d2t
c             enddo
c             y(i,3)=y(i,3)+r2dt8*(x1(i,4)-2.*x1(i,3)+x1(i,2))
c             y(i,kl2)=y(i,kl2)+r2dt8*(x1(i,kles)-2.*x1(i,kl2)
c     1                               +x1(i,kl2-1))/d2t
C
c              IF (ISMG .EQ. 1) THEN
c                 do k=4,kl2-1
c                    Q1_D_V(I,K)=Q1_D_V(I,K)+r2dt2*(Y1(K-1)-Y1(K))
c                 enddo
c                 Q1_D_V(I,3)=Q1_D_V(I,3)
c     1                      +r2dt8*(X1(I,4)-2.*X1(I,3)+X1(I,2))
c                 Q1_D_V(I,KL2)=Q1_D_V(I,KL2)+r2dt8*
c     1                        (X1(I,KLES)-2.*X1(I,KL2)+X1(I,KL2-1))
c              ENDIF
c
c              IF (ISMG .EQ. 2) THEN
c                 do k=4,kl2-1
c                    Q2_D_V(I,K)=Q2_D_V(I,K)+r2dt2*(Y1(K-1)-Y1(K))
c                 enddo
c                 Q2_D_V(I,3)=Q2_D_V(I,3)
c     1                      +r2dt8*(X1(I,4)-2.*X1(I,3)+X1(I,2))
c                 Q2_D_V(I,KL2)=Q2_D_V(I,KL2)+r2dt8*
c     1                        (X1(I,KLES)-2.*X1(I,KL2)+X1(I,KL2-1))
c              ENDIF
c
c          enddo
c
c        ELSE
c 
c               r2dt8=0.001/d2t
c
c          y1(2)=0.
c          y1(kmax)=0.
c          do i=2,iles
c
c            do k=3,kles
c               y1(k)=x1(i,k)-x1(i,k-1)
c            enddo
c            do k=2,kles
c               y(i,k)=y(i,k)-r2dt8*(y1(k)-y1(k+1))
c            enddo
c
c
C
C
c            IF (ISMG .EQ. 1) THEN
c               do k=2,kles
c                 Q1_D_V(I,K)=Q1_D_V(I,K)-0.008*(Y1(K)-Y1(K+1))
c               enddo
c            ENDIF
c            IF (ISMG .EQ. 2) THEN
c               do k=2,kles
c                 Q2_D_V(I,K)=Q2_D_V(I,K)-0.008*(Y1(K)-Y1(K+1))
c               enddo
c            ENDIF

c         enddo
C
C

c        ENDIF
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc
      do 300 k=2,kles
      do 300 i=2,iles
       y1(i)=x(i,k)+eps*(-2.*x(i,k)+x1(i,k))
       x(i,k)=x1(i,k)+y(i,k)*d2t
  300  x1(i,k)=y1(i)
c     ****************************************************************

      ELSE

c     ****************************************************************
       do 401 i=2,iles
  401  y1(i)=Y3(I)
      do 4000 k=3,kmax
        km=k-1
c        A0K=awks*AM1(K)*BSK2(K)
c        A1K=A0K*RDZ
c        A2K=AM(KM)*RRHO(KM)*RD2Z

c       if(iss.eq.1) then
c       do 420 i=2,iles
c        y2(i)=rho1(k)*(-(am1(k)*(x1(i,k)-x1(i,km))
c     1        *(ak(i,k)*(1.+xxw(i,k))+ak(i,km)*(1.+xxw(i,km)))
c     2        +awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,km)))*rdz)
c  420  continue
c       else
       do 440 i=2,iles
       y2(i)=rho1(k)*(-(am1(k)*(x1(i,k)-x1(i,km)+aa1(k)-aa1(km))
     1              *(ak(i,k)*(1.+xxw(i,k))+ak(i,km)*(1.+xxw(i,km)))
     2          +awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,km)))*rdz)
  440  continue
c       endif

       if(isr.eq.0) then
       do 430 i=2,iles
         if (ismg .ne. 4) then
           y17(i)=dlt1*qcl(i,k)+dlt2*qci(i,k)
           y18(i)=dlt1*qcl(i,km)+dlt2*qci(i,km)
         else
           y17(i)=dlt1*qcl_buffer(i,k)+dlt2*qci(i,k)
           y18(i)=dlt1*qcl_buffer(i,km)+dlt2*qci(i,km)
         endif
        if (y17(i) .lt. 1.e-5 .or. y18(i) .lt. 1.e-5) go to 430
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        y2(i)=rho1(k)*(dfc(i,k)
     1                 -awks*am1(k)*bsk2(k)*(x1(i,k)-x1(i,km))*rdz)
  430  continue
       endif
       if (km .eq. kles) then
          do i=2,iles
             y2(i)=0.
          enddo
       endif
       do 470 i=2,iles
        y(i,km)=am(km)*(y1(i)-y2(i))*rd2z*rrho(km)
c        if (km .eq. kles) y(i,km)=am(km)*y1(i)*rd2z*rrho(km)
c
c
c         if (ismg .eq. 1) then
c            q1_d_v(i,km)=q1_d_v(i,km)+am(km)*(y1(i)-y2(i))*rd2z*rrho(km)
c            if (km .eq. kles) 
c     1          q1_d_v(i,km)=q1_d_v(i,km)+am(km)*y1(i)*rd2z*rrho(km)
c         endif
c
c         if (ismg .eq. 2) then
c            q2_d_v(i,km)=q2_d_v(i,km)+am(km)*(y1(i)-y2(i))*rd2z*rrho(km)
c            if (km .eq. kles) 
c     1          q2_d_v(i,km)=q2_d_v(i,km)+am(km)*y1(i)*rd2z*rrho(km)
c         endif
c
c
         y1(i)=y2(i)
  470  continue
 4000 continue

c     ******   horizontal advection terms   **************************

      do 5000 k=2,kles
c         a11k=am(k)*rd4z
         y1(1)=a4h(k)*(x1(3,k)-3.*(x1(2,k)-x1(1,k))-x1(il2,k))

         do i=2,il2
            y1(i)=a4h(k)*(x1(i+2,k)-3.*(x1(i+1,k)-x1(i,k))-x1(i-1,k))
         enddo

         do i=2,il2
c
c
c          if (ismg .eq. 1) then
c            q1_d_h(i,k)=q1_d_h(i,k)+y1(i-1)-y1(i)
c          endif
c          if (ismg .eq. 2) then
c            q2_d_h(i,k)=q2_d_h(i,k)+y1(i-1)-y1(i)
c          endif
c
c
            y(i,k)=y(i,k)+awks*(y1(i-1)-y1(i))
         enddo

c         if (ismg .eq. 1) then
c           q1_d_h(iles,k)=q1_d_h(iles,k)+y1(il2)-y1(1)
c         endif
c         if (ismg .eq. 2) then
c           q2_d_h(iles,k)=q2_d_h(iles,k)+y1(il2)-y1(1)
c         endif

         y(iles,k)=y(iles,k)+awks*(y1(il2)-y1(1))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do 510 i=2,imax
            y1(i)=-(ak(i,k)+ak(i-1,k))*(x1(i,k)-x1(i-1,k))*r2dx2
  510    continue

         do 520 i=2,iles
c
c
c          if (ismg .eq. 1) then
c            q1_d_h(i,k)=q1_d_h(i,k)+y1(i)-y1(i+1)
c          endif
c          if (ismg .eq. 2) then
c            q2_d_h(i,k)=q2_d_h(i,k)+y1(i)-y1(i+1)
c          endif
c
c
            y(i,k)=y(i,k)+y1(i)-y1(i+1)
c
            TRAH(I,K)=Y1(I)-Y1(I+1)
c
            y(i,k)=y(i,k)-a0000*(am(k)*y0(i)*
     1                     (wb(k+1)+wb(k))*(x(i,k+1)-x(i,k-1))*rd4z)
     2                   -a_irf*rfa(k)*x1(i,k)
  520    continue

 5000 continue


      do 405 k=2,kles
      do 405 i=2,iles
  405   xx0(i,k)=y(i,k)
cc
c
       if (ismg.eq.5) then
         call fadvt (x,y,vtp)
          do 495 k=2,kles
          do 495 i=2,iles
  495      xx00(i,k)=y(i,k)-xx0(i,k)
       endif
c
cc
       if (id .eq. 1) then
        do 4150 k=2,kles
        do 415 kc=1,7
         do 425 mt=1,4
         do 425 i=2,iles
           y2(i)=.01*(ww1(i,k)+wb(k)*y0(i))
          ibz(i,mt)=0
          ibz(i,1)=1
           if(ics(i,mt).eq.1) ibz(i,mt)=1
           if(kc.eq.4) go to 435
           if(kc.le.3) go to 445
            if (y2(i).gt.rby(kc)) ibz(i,mt)=0
            go to 435   
  445       if (y2(i).lt.rby(kc)) ibz(i,mt)=0
  435     continue
  425     continue
         do 455 mt=1,4
c          a11=0.
          a22=0.
          a33=0.
         do 465 i=2,iles
          if(ibz(i,mt).eq.1) then
             a33=a33+TRAH(I,K)*a20
c             a11=a11+xx0(i,k)*A20
           if (ismg.gt.2) a22=a22+xx00(i,k)*A20
          endif
  465   continue
         if (ismg.eq.1) sntd(k,mt,kc)=sntd(k,mt,kc)+a33
         if (ismg.eq.2) snqd(k,mt,kc)=snqd(k,mt,kc)+a33
         if (ismg.gt.2) snhd(k,mt,kc)=snhd(k,mt,kc)+a33
         if (ismg.gt.2) snhv(k,mt,kc)=snhv(k,mt,kc)+a22
  455   continue
  415  continue
 4150  continue
      endif
ccc
       call advectn (x,x1)
ccc
       do 590 k=2,kles
       do 590 i=2,iles
  590    x(i,k)=x(i,k)+y(i,k)*dt
      endif
c     ***********************************************************


      if(iss .ne. 0) then
      do 80 k=kt1,kt2
       y1(k)=0.
       y2(k)=0.
       y3(k)=rho(k)*dz/am(k)
   80 continue
      do 82 i=2,iles
      do 82 k=kt1,kt2
       y1(k)=y1(k)+max(x(i,k),0.0)
   82  y2(k)=y2(k)+max(-x(i,k),0.0)
      do 84 k=kt1,kt2
       sqa(k)=sqa(k)+y2(k)*dt/d22t*ril2
   84 continue
      a1=0.
      a2=0.
      do 86 k=kt1,kt2
       a1=a1+y1(k)*y3(k)
   86  a2=a2+y2(k)*y3(k)
      a0=0.
        if(a1.ne.0.) a0=(a1-a2)/a1
      do 88 k=kt1,kt2
      do 88 i=1,imax
       x(i,k)=max(x(i,k),0.0)*a0
   88 continue
      endif
c
      d2t=d22t
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advectn (x,x1)

c     ****   compute advection of t, q, etc

      integer nx,nz
      parameter (NX=514,NZ=43)

      real    x(nx,nz),x1(nx,nz)

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg

      real    u1(nx,nz),w1(nx,nz)
      common/bsat1/ u1
      common/badv1/ w1

      real    uu1(nx,nz),ww1(nx,nz)
      common/b2u/ uu1
      common/b2w/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/dumuw/ umd,vmd,wmd

      real    umod(nx,nz),wmod(nx,nz)
      common/bw0/ umod,wmod

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    aconst
      common/bsave/ aconst

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    TQMD(NX,NZ)
c     real    TQMD1(NX,NZ)
      integer i,k
      real    small

      save TQMD
c     save TQMD1

cc    **************************************************************
        do 10 i=1,imax
         dxr(i)=rdx
   10    dxxt(i)=dt*rdx
       do 20 k=2,kles
        dzr(k)=am(k)*rdz/rho(k)
   20   dzzt(k)=dt*am1(k)*rdz
      do 100 k=1,kmax
      do 100 i=1,imax
         umod(i,k)=0.
         wmod(i,k)=0.
        u1(i,k)=.5*(3.*uu1(i,k)-umd(i,k))
        IF (ISMG .NE. 5) THEN
          W1(I,K)=.5*(3.*WW1(I,K)-WMD(I,K))
        ELSE
          W1(I,K)=WW1(I,K)
        ENDIF
  100 continue

          ACONST=0.0
         IF (ISMG .LE. 2) THEN
          CALL TMAXADV (X,SMALL)
           ACONST=ABS(SMALL*1.125)
         ENDIF

      do 150 k=1,kmax
      do 150 i=1,imax
       IF (ISMG .LE. 2) THEN
        X(I,K)=X(I,K)+ACONST
       ENDIF
       x1(i,k)=x(i,k)
  150 continue
      call fadv (x,u1,w1,0)
       call boundy (x,x1)
cc    *********************************************
      call fadvuw (x,x1,u1,w1,umod,wmod)
c      call boundy (umod,wmod)
      call fadv (x,umod,wmod,0)
       call boundy (x,x1)

      IF (ISMG .LE. 2) THEN

       DO 300 K=1,KMAX
       DO 300 I=1,IMAX
        TQMD(I,K)=ACONST
c        TQMD1(I,K)=TQMD(I,K)
c        U1(I,K)=.5*(3.*UU1(I,K)-UMD(I,K))
c        W1(I,K)=.5*(3.*WW1(I,K)-WMD(I,K))
  300 CONTINUE

c      CALL FADV (TQMD,U1,W1,1)

c       do 50 k=2,kles
c        tqmd(1,k)=tqmd(iles,k)
c  50   tqmd(imax,k)=tqmd(2,k)
c      do 60 i=1,imax
c       tqmd(i,1)=tqmd(i,2)
c  60   tqmd(i,kmax)=tqmd(i,kles)

      DO 400 K=1,KMAX
      DO 400 I=1,IMAX
         X(I,K)=X(I,K)-TQMD(I,K)
         X1(I,K)=X1(I,K)-ACONST
  400 CONTINUE
      ENDIF
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadv (x,u,w,itq)

CC    ****   COMPUTE ADVECTION OF DIFFERENT TRACERS
      integer nx,nz,itq
      PARAMETER (NX=514,NZ=43,ITT=244,NT=2880)

      real x(nx,nz),u(nx,nz),w(nx,nz)

      common/timestat/ ndt_stat,itime_ave,mmave

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      common/bxz/ imax,iles,il2,kmax,kles,kl2

      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg
      common/bcorr/ id

      real    qcl(nx,nz),qrn(nx,nz),qci(nx,nz),qcs(nx,nz),qcg(nx,nz)
      real    qcl1(nx,nz),qrn1(nx,nz),qci1(nx,nz),qcs1(nx,nz),
     1        qcg1(nx,nz),ww1(nx,nz)
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1i/ qci
      common/b1s/ qcs
      common/b1g/ qcg
      common/b2c/ qcl1
      common/b2r/ qrn1
      common/b2i/ qci1
      common/b2s/ qcs1
      common/b2g/ qcg1
      common/b2w/ ww1

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),
     $   y8(nx),y16(nx),y17(nx),y18(nx),x1(nx),x2(nx),x3(nx),x4(nx),
     $   x5(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,y8,y16,y17,y18,x1,x2,x3,x4,x5

      common/bls/ y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/bch/ it(nx),iv(nt),ics(nx,4),ibz(nx,4)

      real    snth(nz,4,7),sntv(nz,4,7),sntd(nz,4,7),snqh(nz,4,7),
     $ snqv(nz,4,7),snqd(nz,4,7),snhh(nz,4,7),snhv(nz,4,7),snhd(nz,4,7)
      common/bsts2/ snth,sntv,sntd,snqh,snqv,snqd,snhh,snhv,snhd

      common/bch1/ rby(7)

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

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,kc,mt
      real    a0,a11,a22

      save
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC    ******   VERTICAL ADVECTION TERMS   ******************************

      a20=ndt_stat*ril2

      DO 1 I=2,ILES
    1  X4(I)=X(I,1)
      DO 1000 K=2,KLES
       A0=AM(K)*RDZ*RRHO(K)
      DO 10 I=2,IMAX
        X1(I)=U(I,K)
        X2(I)=W(I,K)
        X3(I)=W(I,K+1)
        X5(I)=X(I,K)
   10 CONTINUE
      DO 15 I=2,ILES
        Y1(I)=X(I+1,K)
        Y2(I)=X(I,K)
         IF (X1(I+1) .GE. 0.0) Y1(I)=X(I,K)
         IF (X1(I) .GT. 0.0) Y2(I)=X(I-1,K)
        Y3(I)=X(I,K+1)
        Y4(I)=X(I,K)
         IF (X3(I) .GE. 0.0) Y3(I)=X(I,K)
         IF (X2(I) .GT. 0.0) Y4(I)=X4(I)
        Y7(I)=-A0*(RHO1(K+1)*X3(I)*Y3(I)-RHO1(K)*X2(I)*Y4(I))
        Y8(I)=-RDX*(X1(I+1)*Y1(I)-X1(I)*Y2(I))
   15  CONTINUE

      DO I=2,ILES
c
c
c         if (ismg .eq. 1) then
c            q1_g_v(i,k)=q1_g_v(i,k)+y7(i)
c            q1_g_h(i,k)=q1_g_h(i,k)+y8(i)
c         endif
c         if (ismg .eq. 2) then
c            q2_g_v(i,k)=q2_g_v(i,k)+y7(i)
c            q2_g_h(i,k)=q2_g_h(i,k)+y8(i)
c         endif
c
c
         X(I,K)=X(I,K)+(Y7(I)+Y8(I))*DT
         X4(I)=X5(I)
         Y2(I)=.01*(W(I,K)+wb(k)*y0(i))
      ENDDO

      IF (ITQ .EQ. 1) THEN
       DO 20 I=2,ILES
         Y7(I)=-Y7(I)
         Y8(I)=-Y8(I)
   20  CONTINUE
      ENDIF
C     ******   STAT FOR W(QC+QR+QI+QS+QG)             ***************
c
       if(id.eq.1 .or. idq .eq. 1) then
        do 415 kc=1,7
         do 425 mt=1,4
         do 425 i=2,iles
            y16(i)=.01*(ww1(i,k)+wb(k)*y0(i))
          ibz(i,mt)=0
          ibz(i,1)=1
           if(ics(i,mt).eq.1) ibz(i,mt)=1
           if(kc.eq.4) go to 435
           if(kc.le.3) go to 445
            if (y16(i).gt.rby(kc)) ibz(i,mt)=0
            go to 435
  445       if (y16(i).lt.rby(kc)) ibz(i,mt)=0
  435     continue
  425     continue
c
         do 455 mt=1,4
          a11=0.
          a22=0.
         do 465 i=2,iles
          if(ibz(i,mt).eq.1) then
             a11=a11+y8(i)*a20
             a22=a22+y7(i)*a20
          endif
  465   continue
c
         if (ismg.eq.1) snth(k,mt,kc)=snth(k,mt,kc)+a11
         if (ismg.eq.1) sntv(k,mt,kc)=sntv(k,mt,kc)+a22
         if (ismg.eq.2) snqh(k,mt,kc)=snqh(k,mt,kc)+a11
         if (ismg.eq.2) snqv(k,mt,kc)=snqv(k,mt,kc)+a22
         if (ismg.gt.2) snhh(k,mt,kc)=snhh(k,mt,kc)+a11
         if (ismg.gt.2) snhv(k,mt,kc)=snhv(k,mt,kc)+a22

  455   continue
  415  continue
c
      endif
C     ******   end STAT FOR W(QC+QR+QI+QS+QG)             ***********
c
 1000 CONTINUE
      RETURN
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadvuw (x,x1,u,w,umod,wmod)

CC    ****   COMPUTE ADVECTION OF DIFFERENT TRACERS

      integer nx,nz
      PARAMETER (NX=514,NZ=43)

      real    x(nx,nz),u(nx,nz),w(nx,nz),umod(nx,nz),wmod(nx,nz)
      real    x1(nx,nz)

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      common/bxz/ imax,iles,il2,kmax,kles,kl2

      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,epss,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg

      real    uu1(nx,nz),ww1(nx,nz)
      common/b2u/ uu1
      common/b2w/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/dumuw/ umd,vmd,wmd

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),
     $   t(nx),tp(nx),tm(nx),q(nx),qp(nx),qm(nx),y8(nx),y9(nx),y10(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,t,tp,tm,q,qp,qm,y8,y9,y10

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,im,imm,ip,k,km,kmm,kp
      real    a18,ar1,arm,arp,arr
c      real    a0
      real    dtti,dttk,scc
c      real    d1t25

      save
c
       EPS=1.E-10
       SCC=1.0
       A18=1./8.*DT*SCC
c       D1T25=.25*DT*SCC
       DO 5 K=1,KMAX
       DO 5 I=1,IMAX
        UMOD(I,K)=0.
    5   WMOD(I,K)=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 100 K=2,KLES
         KP=K+1
         KM=K-1
c        A0=DZR(K)
        ARR=RRHO(K)
        AR1=RHO1(K)
        ARP=RHO1(KP)
       DO 10 I=2,ILES
        Y4(I)=(X(I,KP)+X(I-1,KP))
        Y5(I)=(X(I,K)+X(I-1,K))
   10   Y6(I)=(X(I,KM)+X(I-1,KM))
       DO 15 I=2,ILES
         IP=I+1
         IM=I-1
         IMM=I-2
         IF(I.EQ.2) IMM=IL2
        IF (Y5(I) .GE. EPS) THEN
           UMOD(I,K)=SCC*(ABS(U(I,K))-U(I,K)*U(I,K)*DXXT(I))
     1                  *(X(I,K)-X(IM,K))/(X(IM,K)+X(I,K)+EPS)
     2      -A18*AM(K)*RDZ*U(I,K)*ARR*(Y4(I)-Y6(I))/(Y4(I)+Y6(I)+EPS)
     3                *(AR1*(W(I,K)+W(IM,K))+ARP*(W(I,KP)+W(IM,KP)))
     4         -U(I,K)*(X(IMM,K)-X(IM,K)-X(I,K)+X(IP,K))/(3.*(X(IMM,K)
     5                  +X(IM,K)+X(I,K)+X(IP,K)+EPS))
c    6         -D1T25*RHO(K)*U(I,K)*(RDX*(U(IP,K)-U(IM,K))
c    7                        +A0*(ARP*(W(I,KP)+W(IM,KP))
c    8                            -AR1*(W(I,K)+W(IM,K))))
        ENDIF
   15  CONTINUE
       UMOD(IMAX,K)=UMOD(2,K)
       UMOD(1,K)=UMOD(ILES,K)
  100 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 300 K=3,KLES
         KP=K+1
         KM=K-1
         KMM=K-2
c       A0=AM1(K)*RDZ*RRHO1(K)
       ARR=RRHO1(K)
       AR1=RHO(K)
       ARM=RHO(KM)
      DO 30 I=2,ILES
        Y4(I)=X(I,K)
        Y5(I)=X(I,K)+X(I,KM)
   30   Y6(I)=X(I,KM)
      DO 35 I=2,ILES
        IP=I+1
        IM=I-1
       IF (Y5(I) .GE. EPS) THEN
         WMOD(I,K)=SCC*(ABS(W(I,K))-W(I,K)*W(I,K)*DZZT(K))
     1                *(X(I,K)-X(I,KM))/(X(I,K)+X(I,KM)+EPS)
     2             -A18*RDX*ARR*W(I,K)
     3                 *(AR1*(U(I,K)+U(IP,K))+ARM*(U(I,KM)+U(IP,KM)))
     4  *(Y4(IP)-Y4(IM)+Y6(IP)-Y6(IM))/(Y4(IP)+Y4(IM)+Y6(IP)+Y6(IM)+EPS)
     5          -W(I,K)*(X(I,KMM)-X(I,KM)-X(I,K)+X(I,KP))/(3.*(X(I,KMM)
     6                  +X(I,KM)+X(I,K)+X(I,KP)+EPS))
c    7       -D1T25*RHO1(K)*W(I,K)*(RDX*ARR*(AR1*(U(IP,K)-U(I,K))
c    8                                       +ARM*(U(IP,KM)-U(I,KM)))
c    9                 +A0*(RHO1(KP)*W(I,KP)-RHO1(KM)*W(I,KM)))
       ENDIF
   35 CONTINUE
       WMOD(IMAX,K)=WMOD(2,K)
       WMOD(1,K)=WMOD(ILES,K)
  300 CONTINUE
      DO 900 I=1,IMAX
       UMOD(I,1)=UMOD(I,2)
       UMOD(I,KMAX)=UMOD(I,KLES)
       WMOD(I,1)=0.0
  900  WMOD(I,KMAX)=0.0
CCCCC   NON-OSCILLATORY OPTION (SMOLARKIEWICZ AND GRABOWSKI, 1990)
      DO 301 K=1,KMAX
      DO 301 I=1,IMAX
       U(I,K)=0.
  301  W(I,K)=0.
      DO 200 K=2,KLES
         KP=K+1
         KM=K-1
          DTTK=DT*AM(K)*RDZ
       DO 20 I=2,ILES
   20   Y5(I)=(X(I,K)+X(I-1,K))
       DO 22 I=2,ILES
         IP=I+1
         IM=I-1
         Y1(I)=0.
         Y2(I)=0.
        IF (Y5(I) .GE. EPS) THEN
         Y1(I)=MAX(X(IM,K),X(I,K),X(IP,K),X1(IM,K),X1(I,K),X1(IP,K))
         Y2(I)=MIN(X(IM,K),X(I,K),X(IP,K),X1(IM,K),X1(I,K),X1(IP,K))
        ENDIF
   22  CONTINUE
       DO 24 I=2,ILES
         IP=I+1
         IM=I-1
         Y3(I)=0.
         Y4(I)=0.
        IF (Y5(I) .GE. EPS) THEN
          DTTI=DT*RDX
        Y3(I)=(Y1(I)-X(I,K))/(DTTI*(MAX(UMOD(I,K),0.)*X(IM,K)
     1                             -MIN(UMOD(IP,K),0.)*X(IP,K))
     2                       +DTTK*(MAX(WMOD(I,K),0.)*X(I,KM)
     3                             -MIN(WMOD(I,KP),0.)*X(I,KP))+EPS)
        Y4(I)=(X(I,K)-Y2(I))/(DTTI*(MAX(UMOD(IP,K),0.)*X(I,K)
     1                             -MIN(UMOD(I,K),0.)*X(I,K))
     2                       +DTTK*(MAX(WMOD(I,KP),0.)*X(I,K)
     3                             -MIN(WMOD(I,K),0.)*X(I,K))+EPS)
        ENDIF
   24  CONTINUE
        Y3(1)=Y3(ILES)
        Y4(1)=Y4(ILES)
       DO 26 I=2,ILES
         IM=I-1
        IF (Y5(I) .GE. EPS) THEN
           U(I,K)=MIN(1.,Y4(IM),Y3(I))*MAX(0.,UMOD(I,K))
     1           +MIN(1.,Y4(I),Y3(IM))*MIN(0.,UMOD(I,K))
        ENDIF
   26  CONTINUE
       U(IMAX,K)=U(2,K)
       U(1,K)=U(ILES,K)
  200 CONTINUE
      DO 28 I=1,NX
       U(I,1)=U(I,2)
   28  U(I,KMAX)=U(I,KLES)
       DO 400 I=2,ILES
         IM=I-1
         IP=I+1
          DTTI=DT*RDX
      DO 40 K=2,KLES
   40   Y5(K)=X(I,K)+X(I,K-1)
       DO 42 K=2,KLES
         KP=K+1
         KM=K-1
         Y1(K)=0.
         Y2(K)=0.
        IF (Y5(K) .GE. EPS) THEN
         Y1(K)=MAX(X(I,KM),X(I,K),X(I,KP),X1(I,KM),X1(I,K),X1(I,KP))
         Y2(K)=MIN(X(I,KM),X(I,K),X(I,KP),X1(I,KM),X1(I,K),X1(I,KP))
        ENDIF
   42  CONTINUE
       DO 44 K=2,KLES
         KP=K+1
         KM=K-1
         Y3(K)=0.
         Y4(K)=0.
        IF (Y5(K) .GE. EPS) THEN
          DTTK=DT*AM(K)*RDZ
         Y3(K)=(Y1(K)-X(I,K))/(DTTI*(MAX(UMOD(I,K),0.)*X(IM,K)
     1                              -MIN(UMOD(IP,K),0.)*X(IP,K))
     2                        +DTTK*(MAX(WMOD(I,K),0.)*X(I,KM)
     3                              -MIN(WMOD(I,KP),0.)*X(I,KP))+EPS)
         Y4(K)=(X(I,K)-Y2(K))/(DTTI*(MAX(UMOD(IP,K),0.)*X(I,K)
     1                              -MIN(UMOD(I,K),0.)*X(I,K))
     2                        +DTTK*(MAX(WMOD(I,KP),0.)*X(I,K)
     3                              -MIN(WMOD(I,K),0.)*X(I,K))+EPS)
        ENDIF
   44  CONTINUE
c       Y3(2)=Y3(3)
c       Y4(2)=Y4(3)
       DO 46 K=3,KLES
         KM=K-1
        IF (Y5(K) .GE. EPS) THEN
           W(I,K)=MIN(1.,Y4(KM),Y3(K))*MAX(0.,WMOD(I,K))
     1           +MIN(1.,Y4(K),Y3(KM))*MIN(0.,WMOD(I,K))
        ENDIF
   46  CONTINUE
  400 CONTINUE
      DO 48 K=3,KLES
       W(IMAX,K)=W(2,K)
       W(1,K)=W(ILES,K)
   48 CONTINUE

      DO 800 K=1,KMAX
      DO 800 I=1,IMAX
       UMOD(I,K)=U(I,K)
  800  WMOD(I,K)=W(I,K)

      RETURN
      END



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadvt (x,y,w)

      implicit none
CC    ****   COMPUTE ADVECTION OF TERMINAL VELOCITY

      integer nx,nz
      PARAMETER (NX=514,NZ=43)

      real    x(nx,nz),y(nx,nz),w(nx,nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),
     $   y8(nx),y16(nx),y17(nx),y18(nx),x1(nx),x2(nx),x3(nx),x4(nx),
     $   x5(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,y8,y16,y17,y18,x1,x2,x3,x4,x5

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a0
      save
      DO 100 K=2,KLES
       A0=AM(K)*RDZ*RRHO(K)
      DO 10 I=2,ILES
        X1(I)=-W(I,K)
   10   X2(I)=-W(I,K+1)
      DO 15 I=2,ILES
        Y1(I)=X(I,K+1)
        Y2(I)=X(I,K)
         IF (X2(I) .GE. 0.0) Y1(I)=X(I,K)
         IF (X1(I) .GT. 0.0) Y2(I)=X(I,K-1)
        Y7(I)=-A0*(RHO1(K+1)*X2(I)*Y1(I)-RHO1(K)*X1(I)*Y2(I))
   15  CONTINUE
      DO 100 I=2,ILES
       Y(I,K)=Y(I,K)+Y7(I)
  100 CONTINUE
      RETURN
      END

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine akcoef (pi)
      parameter (NX=514,NZ=43,ITT=244)
      parameter (nb=nx*nz-6*nx)

      real    pi(nz)

      common/cpsbm/ icps_bm,iexplicit

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
    
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      common/bcor/ irf,iadvh,irfg,idq,ismg

      common/bcor1/ fcor

      common/add_f/ i_four

      common/bstart/ n,isec,nran

      real    dpt1(nx,nz),dqv1(nx,nz),qcl1(nx,nz),qci1(nx,nz),
     $        uuT (nx,nz),wwT (nx,nz),ak1 (nx,nz),ak  (nx,nz)
      common/b2t/ dpt1
      common/b2q/ dqv1
      common/b2c/ qcl1
      common/b2i/ qci1
      common/b1a/ ak
      common/b2u/ uuT
      common/b2w/ wwT
      common/b2a/ ak1
   
      common/badv/ u1(nx,nz)

      common/badv1/ t0(nx),q0(nx),a0(nx),a(nx),a1(nx),a2(nx),ddd(nb)
   
      common/bsat/ uu1(nx,nz)
    
      common/bsat1/ ww1(nx,nz)
    
      common/dumuw/ umd(nx,nz),vmd(nx,nz),wmd(nx,nz)

      common/o4/ a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
   
      common/b4/ tbsk(nz),bskt(nz),bsk2t(nz),bsk(nz),bsk4(nz),
     1           bsit(nz),bsi2t(nz),bsi(nz),bsi4(nz)
   
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
    
      common/ba/ y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),t(nx),
     $        tp(nx),tm(nx),q(nx),qp(nx),qm(nx),y8(nx),y9(nx),y10(nx)
   
      common/bls/ y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)

      common/damp/ rfa(nz),rfa1(nz),cntc(nz),cgwd

c-----------------------------------------------------------------------

      save

c     ********************************
       a0000=0.
         if (IMLIFTING .eq. 1) a0000=1.
       a_irf=0.
         if (irf .eq. 1) a_irf=1.

      alsq=.622*al*al
      cpal=.622*1.61*cp*al
      alcp=al/cp
      cpra=cp*ra
      do 10 k=2,kles
       y1(k)=-am(k)*c1(k)*rd2z
   10  y3(k)=1./ta1(k)
      do 100 k=1,kmax
       y2(k)=1./pi(k)
      do 100 i=1,imax
  100  ak(i,k)=.5*ak(i,k)
c
       IF (IJKADV .EQ. 0) THEN
         do k=1,kmax
         do i=1,imax
            UU1(I,K)=UUT(I,K)
            WW1(I,K)=WWT(I,K)
         enddo
         enddo
       ELSE
         do k=1,kmax
         do i=1,imax
            UU1(I,K)=.5*(3.*UUT(I,K)-UMD(I,K))
            WW1(I,K)=.5*(3.*WWT(I,K)-WMD(I,K))
         enddo
         enddo
       ENDIF
c
c     ********************************
      do 210 i=2,iles
       tm(i)=ta1(1)+dpt1(i,1)
       t(i)=ta1(2)+dpt1(i,2)
       qm(i)=qa1(1)+dqv1(i,1)
       q(i)=qa1(2)+dqv1(i,2)
  210 continue
      do 2000 k=2,kles
        kp=k+1
        km=k-1
c        y13k=y1(k)*y3(k)
       do 250 i=2,iles
         y8(i)=qcl1(i,k)+qci1(i,k)
        tp(i)=ta1(kp)+dpt1(i,kp)
        qp(i)=qa1(kp)+dqv1(i,kp)
        if (y8(i).ge.1.e-5) then
          t0(i)=t(i)*pi(k)
          q0(i)=q(i)
          a0(i)=cpra*t0(i)*t0(i)
          a(i)=(a0(i)+cpal*t0(i)*q0(i))/((a0(i)+alsq*q0(i))*ta1(k))
          a1(i)=tp(i)-tm(i)+alcp*(qp(i)*y2(kp)-qm(i)*y2(km))
          a2(i)=qp(i)-qm(i)+qcl1(i,kp)-qcl1(i,km)+qci1(i,kp)-qci1(i,km)
         u1(i,k)=y1(k)*(a(i)*a1(i)-a2(i))
        else
         u1(i,k)=y1(k)*(tp(i)*(1.+.61*qp(i))-tm(i)*(1.+.61*qm(i)))*y3(k)
        endif
         tm(i)=t(i)
         t(i)=tp(i)
         qm(i)=q(i)
         q(i)=qp(i)
  250  continue
 2000 continue

      do 310 i=1,imax
       qm(i)=ak1(i,1)*ak1(i,1)
       tm(i)=ak1(i,2)*ak1(i,2)
  310 continue

      do 3000 k=2,kles
        kp=k+1
        km=k-1
        a11k=am(k)*rd4z
c        c11k=2.*coef(k)*am(k)*am(k)*rdz2
       do 320 i=1,imax
  320   tp(i)=ak1(i,kp)*ak1(i,kp)
       do 330 i=2,iles
         y9(i)=a0000*y0(i)*wb(k)
         y10(i)=a0000*y0(i)*wb(kp)
       y1(i)=-c2(k)*tm(i)+coef(k)*
     1   (((ww1(i+1,kp)+ww1(i+1,k)-ww1(i-1,kp)-ww1(i-1,k))*rd4x
     2     +(uu1(i+1,kp)+uu1(i,kp)-uu1(i+1,km)-uu1(i,km))*a11k)**2
     3    +2.*(uu1(i+1,k)-uu1(i,k))**2*rdx2)
     4   +(tm(i+1)+tm(i-1)-2.*tm(i))*r2dx2
     5   +am(k)*(am1(kp)*(tp(i)-tm(i))-am1(k)*(tm(i)-qm(i)))*r2dz2
     6   +2.*coef(k)*am(k)*am(k)*rdz2
     7                 *(ww1(i,kp)-ww1(i,k)+y10(i)-y9(i))**2
       u1(i,k)=u1(i,k)+y1(i)
       u1(i,k)=u1(i,k)-a0000*(a11k*y0(i)*
     1                   (wb(k+1)+wb(k))*(ak(i,k+1)-ak(i,k-1)))
  330  continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       do 390 i=1,imax
        qm(i)=tm(i)
  390   tm(i)=tp(i)
 3000 continue
c
      if (ijkadv .eq. 1) then
c
       awks=.5
c
      do 5000 k=2,kles
        kp=k+1
        km=k-1
c         y9(i)=y0(i)*wb(k)
c         y10(i)=y0(i)*wb(kp)
c          IF (IMLIFTING .EQ. 0) THEN
           do i=2,iles
             y9(i)=0.0
             y10(i)=0.0
           enddo
c          ENDIF
        a22z=AM(K)*RRHO(K)*RD2Z
        ammkp=awks*am1(kp)*rho1(kp)*bsk(kp)
        ammk=awks*am1(k)*rho1(k)*bsk(k)
       do 520 i=2,iles
        u1(i,k)=u1(i,k)-a22z*(
     1   rho1(kp)*y10(i)*(ak(i,k)+ak(i,kp))
     2   -rho1(k)*y9(i)*(ak(i,k)+ak(i,km))
     3   -2.*(ammkp*(ak1(i,kp)-ak1(i,k))-ammk*(ak1(i,k)-ak1(i,km)))*rdz)
        U1(I,K)=U1(I,K)-a_irf*RFA(K)*AK1(I,K)
  520  continue
        y2(1)=a4k(k)*(ak1(3,k)-3.*(ak1(2,k)-ak1(1,k))-ak1(il2,k))
        do 540 i=2,il2
         y2(i)=a4k(k)*(ak1(i+2,k)-3.*(ak1(i+1,k)-ak1(i,k))-ak1(i-1,k))
  540   continue
        do 570 i=2,il2
  570    u1(i,k)=u1(i,k)+awks*(y2(i-1)-y2(i))
        u1(iles,k)=u1(iles,k)+awks*(y2(il2)-y2(1))

        do i=2,imax
           y1(i)=-(ak(i,k)+ak(i-1,k))*(ak1(i,k)-ak1(i-1,k))*r2dx2
        enddo
        do i=2,iles
           u1(i,k)=u1(i,k)+y1(i)-y1(i+1)
        enddo

 5000 continue
c
       do i=2,iles
          y1(i)=0.
       enddo
      do 4000 k=3,kmax
        km=k-1
       do i=2,iles
          y2(i)=rho1(k)*(-(am1(k)*(ak1(i,k)-ak1(i,km))
     1                  *(ak(i,k)+ak(i,km)))*rdz)
       enddo
       if (km .eq. kles) then
          do i=2,iles
             y2(i)=0.
          enddo
       endif
       do i=2,iles
        u1(i,km)=u1(i,km)+am(km)*(y1(i)-y2(i))*rd2z*rrho(km)
         y1(i)=y2(i)
       enddo
 4000 continue
cc
       call advectak (ak,ak1)
cc
      else
cc
      DO 6000 K=2,KLES
        KP=K+1
        KM=K-1
c        a22z=AM(K)*RRHO(K)*RD2Z
        ammkp=AM1(KP)*RHO1(KP)*BSK(KP)
        ammk=AM1(K)*RHO1(K)*BSK(K)
c         y9(i)=y0(i)*wb(k)
c         y10(i)=y0(i)*wb(kp)
c        IF (IMLIFTING .EQ. 0) THEN
           do i=2,iles
             y9(i)=0.0
             y10(i)=0.0
           enddo
c        ENDIF
       DO 610 I=2,ILES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          u1(i,k)=u1(i,k)-am(k)*rrho(k)*rd2z*(
     1      rho1(kp)*(ww1(i,kp)+y10(i))*(ak(i,k)+ak(i,kp))-2.
     2      *(ammkp*(ak1(i,kp)-ak1(i,k))-ammk*(ak1(i,k)-ak1(i,km)))*rdz        
     3      -rho1(k)*(ww1(i,k)+y9(i))*(ak(i,k)+ak(i,km)))
          U1(I,K)=U1(I,K)-a_irf*RFA(K)*AK1(I,K)
  610  CONTINUE
       if(iadvh .eq. 2) then
        do 620 i=2,iles
         u1(i,k)=u1(i,k)
     1   -(uu1(i+1,k)*(ak(i+1,k)+ak(i,k))-uu1(i,k)*(ak(i,k)+ak(i-1,k)))
c    2   *rd2x+bsi(k)*rdx2*(ak1(i+1,k)-2.*ak1(i,k)+ak1(i-1,k))
     2   *rd2x
  620   continue
        y2(1)=a4k(k)*(ak1(3,k)-3.*(ak1(2,k)-ak1(1,k))-ak1(il2,k))
        do 615 i=2,il2
         y2(i)=a4k(k)*(ak1(i+2,k)-3.*(ak1(i+1,k)-ak1(i,k))-ak1(i-1,k))
  615   continue
        do 625 i=2,il2
  625    u1(i,k)=u1(i,k)+y2(i-1)-y2(i)
        u1(iles,k)=u1(iles,k)+y2(il2)-y2(1)
       else
c     ***   4-th order horizontal advection scheme   ************
        y2(1)=uu1(2,k)*(d58x*(ak(2,k)+ak(1,k))-d16x*(ak(3,k)+ak(il2,k)))
     1        +a4k(k)*(ak1(3,k)-3.*(ak1(2,k)-ak1(1,k))-ak1(il2,k))
        do 630 i=2,il2
         y2(i)=uu1(i+1,k)*(d58x*(ak(i+1,k)+ak(i,k))
     1                     -d16x*(ak(i+2,k)+ak(i-1,k)))
     2         +a4k(k)*(ak1(i+2,k)-3.*(ak1(i+1,k)-ak1(i,k))-ak1(i-1,k))
  630   continue
        y3(1)=-d48x*uu1(1,k)*(ak(1,k)+ak(il2,k))
        do 640 i=2,imax
  640    y3(i)=-d48x*uu1(i,k)*(ak(i,k)+ak(i-1,k))
        do 650 i=2,il2
  650    u1(i,k)=u1(i,k)+y2(i-1)-y2(i)+y3(i-1)-y3(i+2)
        u1(iles,k)=u1(iles,k)+y2(il2)-y2(1)+y3(il2)-y3(3)
       endif
 6000 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-13-97)

c        IF (I_FOUR .EQ. 1) THEN
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c           do i=2,iles
c             do k=3,kl2-1
c               y1(k)=(ak1(i,k+2)-3.*(ak1(i,k+1)-ak1(i,k))-ak1(i,k-1))
c             enddo
c             do k=4,kl2-1
c                u1(i,k)=u1(i,k)+r2dt2*(y1(k-1)-y1(k))/d2t
c             enddo
c             u1(i,3)=u1(i,3)+r2dt8*(ak1(i,4)-2.*ak1(i,3)+ak1(i,2))
c             u1(i,kl2)=u1(i,kl2)+r2dt8*(ak1(i,kles)-2.*ak1(i,kl2)
c     1                                 +ak1(i,kl2-1))
c          enddo
c
c
c        ELSE
c 
c          y1(2)=0.
c          y1(kmax)=0.

c               r2dt8=0.008/d2t

c            do i=2,iles

c            do k=3,kles
c               y1(k)=ak1(i,k)-ak1(i,k-1)
c            enddo
c            do k=2,kles
c               u1(i,k)=u1(i,k)-r2dt8*(y1(k)-y1(k+1))
c            enddo

c          enddo
C
C
cc        ENDIF
c
        IF (I_FOUR .EQ. 1) THEN
           R2DT2=0.005/DT
           DO I=2,ILES
             DO K=2,KLES
               Y1(K)=AK1(I,K+1)-2.*AK1(I,K)+AK1(I,K-1)
             ENDDO
             Y1(1)=Y1(2)
             Y1(KMAX)=Y1(KLES)
             DO K=2,KLES
                U1(I,K)=U1(I,K)-R2DT2*(Y1(K+1)-2.*Y1(K)+Y1(K-1))
             ENDDO
           ENDDO
        ENDIF
cc
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      endif
c
c     ***********************************************************
      do 400 k=2,kles
       a11=tbsk(k)
      do 400 i=2,iles
       IF (IJKADV .EQ. 0) THEN
         Y1(I)=AK(I,K)+EPS*(-2.*AK(I,K)+AK1(I,K))
        AK(I,K)=AK1(I,K)+U1(I,K)*D2T
        AK1(I,K)=Y1(I)+EPS*AK(I,K)
       ELSE
        AK(I,K)=AK(I,K)+U1(I,K)*DT
c        ak1(i,k)=ak(i,k)
       ENDIF
        ak(i,k)=min(a11, max(ak(i,k), 0.0))
        ak1(i,k)=min(a11, max(ak1(i,k), 0.0))
  400 continue
      call boundy (ak,ak1)
      return
      end


Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine akcoefRH 
      integer nx,nz
      parameter (NX=514,NZ=43)

      integer LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    uuT (nx,nz),wwT(nx,nz),ak1(nx,nz),ak(nx,nz)

      common/b1a/ ak
      common/b2u/ uuT
      common/b2w/ wwT
      common/b2a/ ak1

      real    u1(nx,nz)
      common/badv/ u1

      real    uu1(nx,nz)
      common/bsat/ uu1

      real    ww1(nx,nz)
      common/bsat1/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/dumuw/ umd,vmd,wmd

      real    tbsk(nz),bskt(nz),bsk2t(nz),bsk(nz),bsk4(nz),
     1        bsit(nz),bsi2t(nz),bsi(nz),bsi4(nz)
      common/b4/ tbsk,bskt,bsk2t,bsk,bsk4,bsit,bsi2t,bsi,bsi4

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef22(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef22,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),t(nx),
     $        tp(nx),tm(nx),q(nx),qp(nx),qm(nx),y8(nx),y9(nx),y10(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,t,tp,tm,q,qp,qm,y8,y9,y10

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,km,kp
      real    a111,coef(nz)

      save

c     ********************************
      DO K=1,KMAX
         Y1(K)=DZ/AM(K)
        COEF(K)=.5*CK*CK*Y1(K)*Y1(K)
      enddo
      do 100 k=1,kmax
      do 100 i=1,imax
       IF (IJKADV .EQ. 0) THEN
         UU1(I,K)=UUT(I,K)
         WW1(I,K)=WWT(I,K)
       ELSE
         UU1(I,K)=.5*(3.*UUT(I,K)-UMD(I,K))
         WW1(I,K)=.5*(3.*WWT(I,K)-WMD(I,K))
       ENDIF
  100  ak(i,k)=.5*ak(i,k)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ***   CONVENTIONAL MIXING COEFFICIENT   **************************
c      do 200 k=2,kles
c        kp=k+1
c        km=k-1
c        a111=am(k)*rd2z*rrho(k)*(rho1(kp)-rho1(k))
c        do i=2,iles
c          y1(i)=(ww1(i+1,kp)+ww1(i+1,k)-ww1(i-1,kp)-ww1(i-1,k))*rd4x
c     1         +(uu1(i+1,kp)+uu1(i,kp)-uu1(i+1,km)-uu1(i,km))*rd4z*am(k)
c     2         -a111*(ww1(i,kp)+ww1(i,k))
c          u1(i,k)=coef(k)*y1(i)
c       enddo
c  200  continue
c      do 400 k=2,kles
c       a111=tbsk(k)
c       do i=2,iles
c         ak(i,k)=abs( u1(i,k) )
c          if (ak(i,k) .ge. a111) ak(i,k)=a111
c         ak1(i,k)=ak(i,k)
c       enddo
c  400  continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     ***   SCHLESINGER'S MIXING COEFFICIENT   *************************
       A13=1./3.
      do 200 k=2,kles
         kp=k+1
         km=k-1
        a111=a13*am(k)*rd2z*rrho(k)*(rho1(kp)-rho1(k))
        do i=2,iles
           tp(i)=a111*(ww1(i,kp)+ww1(i,k))
          y1(i)=(ww1(i+1,kp)+ww1(i+1,k)-ww1(i-1,kp)-ww1(i-1,k))*rd4x
     1         +(uu1(i+1,kp)+uu1(i,kp)-uu1(i+1,km)-uu1(i,km))*rd4z*am(k)
          y1(i)=y1(i)*y1(i)

          y2(i)=(uu1(i+1,k)-uu1(i,k))*rdx+tp(i)
           y2(i)=y2(i)*y2(i) 
          y3(i)=am(k)*(ww1(i,kp)-ww1(i,k))*rdz+tp(i)
           y3(i)=y3(i)*y3(i)
          u1(i,k)=coef(k)*(2.*(y2(i)+y3(i))+y1(i))
        enddo
  200 CONTINUE

      do 400 k=2,kles
         a111=tbsk(k)
        do i=2,iles
         ak(i,k)=sqrt(u1(i,k))
          if (ak(i,k) .ge. a111) ak(i,k)=a111
         ak1(i,k)=ak(i,k)
        enddo
  400 continue

      call boundy (ak,ak1)

      return
      end




Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fadvak (x,u,w)
      implicit none
CC    ****   COMPUTE ADVECTION OF DIFFERENT TRACERS
      integer nx,nz
      PARAMETER (NX=514,NZ=43)
      real    x(nx,nz),u(nx,nz),w(nx,nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr
      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1
      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx),
     $   y8(nx),y16(nx),y17(nx),y18(nx),x1(nx),x2(nx),x3(nx),x4(nx),
     $   x5(nx)
      common/ba/ y1,y2,y3,y4,y5,y6,y7,y8,y16,y17,y18,x1,x2,x3,x4,x5
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a0
      save
CC    ******   VERTICAL ADVECTION TERMS   **************************
      DO 1 I=2,IMAX
    1  X4(I)=X(I,1)
      DO 100 K=2,KLES
       A0=AM(K)*RDZ*RRHO(K)
      DO 10 I=2,IMAX
        X1(I)=U(I,K)
        X2(I)=W(I,K)
        X3(I)=W(I,K+1)
   10   X5(I)=X(I,K)
      DO 15 I=2,ILES
        Y1(I)=X(I+1,K)
        Y2(I)=X(I,K)
         IF (X1(I+1) .GE. 0.0) Y1(I)=X(I,K)
         IF (X1(I) .GT. 0.0) Y2(I)=X(I-1,K)
        Y3(I)=X(I,K+1)
        Y4(I)=X(I,K)
         IF (X3(I) .GE. 0.0) Y3(I)=X(I,K)
         IF (X2(I) .GT. 0.0) Y4(I)=X4(I)
        Y7(I)=-A0*(RHO1(K+1)*X3(I)*Y3(I)-RHO1(K)*X2(I)*Y4(I))
        Y8(I)=-RDX*(X1(I+1)*Y1(I)-X1(I)*Y2(I))
   15  CONTINUE
      DO 100 I=2,ILES
       X(I,K)=X(I,K)+(Y7(I)+Y8(I))*DT
       X4(I)=X5(I)
  100 CONTINUE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advwuv (nudge)
      integer nx,nz,nx10,itt
      parameter (NX=514,NZ=43,ITT=244,ITT3=ITT*3)
      parameter (nx10=10*nx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      COMMON/SFLUXS/ SUW,SVW,SWT,SWQ
      integer itoga,ISFC,ICE,ICE2
      common/itoga/ itoga,ISFC,ICE,ICE2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/cpsbm/ icps_bm,iexplicit

      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      integer irf,iadvh,irfg,idq,ismg
      common/bcor/ irf,iadvh,irfg,idq,ismg

      real    fcor
      common/bcor1/ fcor

      common/add_f/ i_four

      integer n,isec,nran
      common/bstart/ n,isec,nran

      real    dpt(nx,nz),dqv(nx,nz),qcl(nx,nz),qrn(nx,nz),
     $        qci(nx,nz),qcs(nx,nz),qcg(nx,nz),
     $        u  (nx,nz),v  (nx,nz),w  (nx,nz),ak (nx,nz),
     $        uu1(nx,nz),vv1(nx,nz),ww1(nx,nz)
      common/b1t/ dpt
      common/b1q/ dqv
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1i/ qci
      common/b1s/ qcs
      common/b1g/ qcg
      common/b1u/ u
      common/b1v/ v
      common/b1w/ w
      common/b1a/ ak
      common/b2u/ uu1
      common/b2v/ vv1
      common/b2w/ ww1

      real    u1(nx,nz)
      common/bsat/ u1

      real    w1(nx,nz)
      common/bsat1/ w1

      real    v1(nx,nz)
      common/badv/ v1

      real    a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
      common/o4/ a4k,a4h,d58x,d16x,d48x,d516x,d32x,d96x

      real    tbsk(nz),bskt(nz),bskt2(nz),bsk(nz),bsk4(nz),
     1        bsit(nz),bsit2(nz),bsi(nz),bsi4(nz)
      common/b4/ tbsk,bskt,bskt2,bsk,bsk4,bsit,bsit2,bsi,bsi4

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/bb6/ tls1,tls2,qls1,qls2,tls3,tls4,qls3,qls4,sft,sfq,wbt,
     $            wb_6h,ub_6h,ubt,q1_6h,q1t,q2_6h,q2t,vb_6h,vbt
      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx10)
      common/ba/ y1,y2,y3,y4,y5,y6,y7
      real    y11(nx),y22(nx),y33(nx),y44(nx),y55(nx),y66(nx)
      real    y0(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/bls/y0,ths,qs,ts,pss
      real    rfa(nz),rfa1(nz),cntc(nz),cgwd
      common/damp/ rfa,rfa1,cntc,cgwd

      common/bls3/ factor_nuding
      common/bls4/ ubi(nz),vbi(nz),ub_2h(nz,itt3),vb_2h(nz,itt3)
C
CC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFF_VX(NX,NZ) :  V-TURBULENT IN X-COMPNT                        C
C     DIFF_VZ(NX,NZ) :  V-TURBULENT IN Z-COMPNT                        C
C     GRID_VX(NX,NZ) :  V-GRID SCALE TRANSPT IN X-COMPNT               C
C     GRID_VZ(NX,NZ) :  V-GRID SCALE TRANSPT IN Z-COMPNT               C
C     DIFF_NV(NX,NZ) :  V-NUMERICAL FILTER                             C
C     V_LARGE(NX,NZ) :  LARGE-SCALE FORCING IN V WIND                  C
C     PRE_V(NX,NZ)   :  HORIZONTAL PREES FORCING                       C
C     DT_VWIND(NX,NZ):  LOCAL TIME CHANGE TERM                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     DIFF_UX(NX,NZ) :  U-TURBULENT IN X-COMPNT                        C
C     DIFF_UZ(NX,NZ) :  U-TURBULENT IN Z-COMPNT                        C
C     GRID_UX(NX,NZ) :  U-GRID SCALE TRANSPT IN X-COMPNT               C
C     GRID_UZ(NX,NZ) :  U-GRID SCALE TRANSPT IN Z-COMPNT               C
C     DIFF_NU(NX,NZ) :  U-NUMERICAL FILTER                             C
C     U_LARGE(NX,NZ) :  LARGE-SCALE FORCING IN U WIND                  C
C     PRE_U(NX,NZ)   :  HORIZONTAL PREES FORCING                       C
C     DT_UWIND(NX,NZ):  LOCAL TIME CHANGE TERM                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real            DIFF_VX(NX,NZ),DIFF_VZ(NX,NZ),GRID_VX(NX,NZ),
     1                GRID_VZ(NX,NZ),DIFF_NV(NX,NZ),V_LARGE(NX,NZ),
     2                PRE_V(NX,NZ),DT_VWIND(NX,NZ)
      COMMON/DEBUG_V/ DIFF_VX,DIFF_VZ,GRID_VX,
     1                GRID_VZ,DIFF_NV,V_LARGE,
     2                PRE_V,DT_VWIND

      real            DIFF_UX(NX,NZ),DIFF_UZ(NX,NZ),GRID_UX(NX,NZ),
     1                GRID_UZ(NX,NZ),DIFF_NU(NX,NZ),U_LARGE(NX,NZ),
     2                PRE_U(NX,NZ),DT_UWIND(NX,NZ)
      COMMON/DEBUG_U/ DIFF_UX,DIFF_UZ,GRID_UX,
     1                GRID_UZ,DIFF_NU,U_LARGE,
     2                PRE_U,DT_UWIND

      real UTEM_DT(NX,NZ),VTEM_DT(NX,NZ)
      COMMON/DEBUG_UUVV/ UTEM_DT,VTEM_DT

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,km,kmm,kp,im
      real    a0000,c33k,rdxz4,ub2,BARFCT,rd2t,ub12,vb12

      save

c     ******   compute w1   ********************************************
      A0000=1.
       IF (IMLIFTING .EQ. 0) A0000=0.
       a_irf=0.
         if (irf .eq. 1) a_irf=1.
c
      do 10 i=2,iles
       y1(i)=0.
       y3(i)=0.
   10 continue
      do 1000 k=2,kles
        kp=k+1
c        a11k=AM(K)*RDZ2
c        a22k=am1(k)*rrho1(k)
       do 100 i=2,iles
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDANCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        y2(i)=rho(k)*(w(i,kp)+w(i,k))**2*rd4z
        y4(i)=rho(k)*(-am(k)*
     1        (2.*ak(i,k)+bsk(k))
c    1        (2.*ak(i,k)+.5*(bsk(k)+bsk(kp)))
     2        *(ww1(i,kp)-ww1(i,k))*rdz2
     3        +c3(k)*ak(i,k)**2*rdz)
       w1(i,k)=am1(k)*(y1(i)-y2(i)+y3(i)-y4(i))*rrho1(k)
        y1(i)=y2(i)
        y3(i)=y4(i)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  100  continue
 1000 continue
c      do 115 i=2,iles
c        w1(i,kmax)=am1(kmax)*(y2(i)+y4(i))*rrho1(kmax)
c  115 continue

cccccc  tao (11-13-97)
        if (ijkadv .eq. 0 .and. I_FOUR .EQ. 1) then
cc        IF (I_FOUR .EQ. 1) THEN
c
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c
c          do i=2,iles
c             do k=3,kl2
c               y1(k)=(ww1(i,k+2)-3.*(ww1(i,k+1)-ww1(i,k))-ww1(i,k-1))
c             enddo
c             do k=4,kl2
c               w1(i,k)=w1(i,k)+r2dt2*(y1(k-1)-y1(k))
c             enddo
c             w1(i,3)=w1(i,3)+r2dt8*(ww1(i,4)-2.*ww1(i,3))
c             w1(i,kles)=w1(i,kles)+r2dt8*(-2.*ww1(i,kles)+ww1(i,kl2))
c          enddo
c
c          IF (I_FOUR .EQ. 1) THEN
c
               r2dt2=0.005/dt
c
          do i=2,iles
             do k=2,kles
               y1(k)=ww1(i,k+1)-2.*ww1(i,k)+ww1(i,k-1)
             enddo
             y1(2)=y1(3)
             y1(kmax)=y1(kles)
             do k=3,kles
               w1(i,k)=w1(i,k)-r2dt2*(y1(k+1)-2.*y1(k)+y1(k-1))
             enddo
          enddo

c
c
c        ELSE
c
c          y1(2)=0.
c          y1(kmax)=0.
c               r2dt8=0.008/d2t
c          do i=2,iles
c
c            do k=3,kles
c               y1(k)=ww1(i,k)-ww1(i,k-1)
c            enddo
c            do k=3,kles
c               w1(i,k)=w1(i,k)-r2dt8*(y1(k)-y1(k+1))
c            enddo

c          enddo
C
C

cc        ENDIF
        endif
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       IF (IMLIFTING .EQ. 1) THEN
        do 30 k=2,kmax
   30    y1(k)=0.
        do 40 i=2,iles
        do 40 k=2,kmax
   40    y1(k)=y1(k)+w1(i,k)
        do 45 k=2,kmax
   45     y1(k)=y1(k)*ril2
        do 50 k=2,kmax
        do 50 i=2,iles
   50    w1(i,k)=w1(i,k)-y1(k)
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        rdxz4=rd4z*rdx
      do 1500 k=3,kles
        kp=k+1
        km=k-1
        kmm=k-2
c        a11k=AM1(K)*RD2Z
c        a22k=am1(k)*rdxz4
       if (iadvh .eq. 2) then
       do 120 i=2,imax
         y6(i)=ak(i,k)+ak(i,km)+ak(i-1,k)+ak(i-1,km)
        y1(i)=(u(i,k)+u(i,km))*(w(i,k)+w(i-1,k))*rd4x
c    1        -(y6(i)+bsi4(k))*(ww1(i,k)-ww1(i-1,k))*r4dx2
     1        -y6(i)*(ww1(i,k)-ww1(i-1,k))*r4dx2
     2        -am1(k)*y6(i)*(uu1(i,k)-uu1(i,km))*rd4z*rdx
  120  continue
       do 130 i=2,iles
  130   w1(i,k)=w1(i,k)+y1(i)-y1(i+1)
c       y1(1)=(ww1(3,k)-3.*(ww1(2,k)-ww1(1,k))-ww1(il2,k))
c       do 125 i=2,il2
c  125   y1(i)=(ww1(i+2,k)-3.*(ww1(i+1,k)-ww1(i,k))-ww1(i-1,k))
c       do 135 i=2,il2
c  135   w1(i,k)=w1(i,k)+.5*(a4k(k)+a4k(km))*(y1(i-1)-y1(i))
c       w1(iles,k)=w1(iles,k)+.5*(a4k(k)+a4k(km))*(y1(il2)-y1(1))
       y1(1)=a4k(k)*(ww1(3,k)-3.*(ww1(2,k)-ww1(1,k))-ww1(il2,k))
       do 125 i=2,il2
  125   y1(i)=a4k(k)*(ww1(i+2,k)-3.*(ww1(i+1,k)-ww1(i,k))-ww1(i-1,k))
       do 135 i=2,il2
  135   w1(i,k)=w1(i,k)+y1(i-1)-y1(i)
       w1(iles,k)=w1(iles,k)+y1(il2)-y1(1)
       else
c     ***   4-th order horizontal advection
       y3(1)=-d96x*(u(1,kmm)+u(1,kp))*(w(il2,k)+w(1,k))
       do 140 i=2,imax
         y6(i)=ak(i,k)+ak(i,km)+ak(i-1,k)+ak(i-1,km)
        y1(i)=-y6(i)*((ww1(i,k)-ww1(i-1,k))*r4dx2
     1        +am1(k)*(uu1(i,k)-uu1(i,km))*rd4z*rdx)
        y3(i)=-d96x*(u(i,kmm)+u(i,kp))*(w(i-1,k)+w(i,k))
  140  continue
       do 150 i=2,iles
  150   w1(i,k)=w1(i,k)+y1(i)-y1(i+1)
       y1(1)=(u(2,km)+u(2,k))*(d516x*(w(1,k)+w(2,k))
     1                         -d32x*(w(il2,k)+w(3,k)))
     2       +a4k(k)*
c    2       +.5*(a4k(k)+a4k(km))*
     3               (ww1(3,k)-3.*(ww1(2,k)-ww1(1,k))-ww1(il2,k))
       do 160 i=2,il2
        y1(i)=(u(i+1,km)+u(i+1,k))*(d516x*(w(i,k)+w(i+1,k))
     1                              -d32x*(w(i-1,k)+w(i+2,k)))
     2        +a4k(k)*
c    2        +.5*(a4k(k)+a4k(km))*
     3             (ww1(i+2,k)-3.*(ww1(i+1,k)-ww1(i,k))-ww1(i-1,k))
  160  continue
       do 170 i=2,il2
  170   w1(i,k)=w1(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)
       w1(iles,k)=w1(iles,k)+y1(il2)-y1(1)+y3(il2)-y3(3)
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO 175 I=2,ILES
  175   W1(I,K)=W1(I,K)
c    1         -a000*a11k*Y0(I)*(WB(KP)*W(I,KP)-WB(KM)*W(I,KM))
     1         -a0000*am1(k)*y0(i)*(wb(kp)*w(i,kp)-wb(km)*w(i,km))*rd2z
 1500 continue
      do i=2,iles
         w1(i,2)=0.
         w1(i,kmax)=0.
      enddo
c     ******   compute buoyancy
      do 180 i=2,iles
        y1(i)=490.*(dpt(i,1)/tb(1)+.61*dqv(i,1)-qcl(i,1)
     1        -qrn(i,1)-qci(i,1)-qcs(i,1)-qcg(i,1)+a0000*sq(1))
  180  continue
      do 1800 k=2,kmax
       do 190 i=2,iles
         y2(i)=490.*(dpt(i,k)/tb(k)+.61*dqv(i,k)-qcl(i,k)
     1         -qrn(i,k)-qci(i,k)-qcs(i,k)-qcg(i,k)+a0000*sq(k))
        w1(i,k)=w1(i,k)+y1(i)+y2(i)
        w1(i,k)=w1(i,k)-a_irf*rfa1(k)*ww1(i,k)
  190    y1(i)=y2(i)
 1800 continue
c     ******   compute u1   *******************************************
      do i=2,iles
         y1(i)=0.
         y2(i)=0.
         y3(i)=0.
         y4(i)=0.
         y5(i)=0.
         y6(i)=0.
         y7(i)=0.
         y11(i)=0.
         y22(i)=0.
         y33(i)=0.
         y44(i)=0.
         y55(i)=0.
         y66(i)=0.
      enddo
C
C-----put the surface momentum fluxes in here
C
      IF (ITOGA .EQ. 1 .AND. ISFC .EQ. 1) THEN
        BARFCT=.5*RHO1(2)*100.*100.*RDZ  !convert from mks to cgs
        DO I=2,ILES
            im=i-1
            if (i.eq.2) im=iles
          y3(i)=(suw(i)+suw(im))*barfct
          y33(i)=(svw(i)+svw(im))*barfct
        ENDDO
      ENDIF
       rdxz4=rd4x*rdz
      do 2000 k=3,kles
       km=k-1
c        a11_4z2=am1(k)*r4dz2
c        a11_rhoz=rho1(k)*rd4z
c        a11z=am(km)*rrho(km)
c        BASE_DIFF=rho1(k)*am1(k)*bsk4(k)*r4dz2
c
        ub12=ub(k)+ub(km)
        vb12=vb(k)+vb(km)
c
       do 200 i=2,iles
         y7(i)=ak(i,k)+ak(i,km)+ak(i-1,k)+ak(i-1,km)
c
c         y2(i)=rho1(k)*((w(i,k)+w(i-1,k))*(u(i,k)+u(i,km))*rd4z
c     1        -am1(k)*r4dz2*(y7(i)*(uu1(i,k)-uu1(i,km))
c     2        +bsk4(k)*(uu1(i,k)-uu1(i,km)-ub1(k)+ub1(km)))
c     3        -y7(i)*(ww1(i,k)-ww1(i-1,k))*rdxz4)
c
c         y22(i)=rho1(k)*((w(i,k)+w(i-1,k))*(v(i,k)+v(i,km))*rd4z
c     1        -am1(k)*r4dz2*(y7(i)*(vv1(i,k)-vv1(i,km))
c     2        +bsk4(k)*(vv1(i,k)-vv1(i,km)-vb1(k)+vb1(km))))
c

         y2(i)=rho1(k)*(w(i,k)+w(i-1,k))*(u(i,k)+u(i,km))*rd4z
         Y6(I)=RHO1(K)*(W(I,K)+W(I-1,K))*(U(I,K)+U(I,KM)-UB12)*RD4Z

c         y4(i)=-rho1(k)*(am1(k)*r4dz2*y7(i)*(uu1(i,k)-uu1(i,km))
c     1                               +y7(i)*(ww1(i,k)-ww1(i-1,k))*rdxz4)
c         y6(i)=-rho1(k)*am1(k)*r4dz2*bsk4(k)*(uu1(i,k)-uu1(i,km)
c     1                                       -ub1(k)+ub1(km))

         y4(i)=rho1(k)*(-am1(k)*r4dz2*(y7(i)*(uu1(i,k)-uu1(i,km))
     1                  +bsk4(k)*(uu1(i,k)-uu1(i,km)-ub1(k)+ub1(km)))
     2                  -y7(i)*(ww1(i,k)-ww1(i-1,k))*rdxz4)

c

         y22(i)=rho1(k)*(w(i,k)+w(i-1,k))*(v(i,k)+v(i,km))*rd4z
         Y66(I)=RHO1(K)*(W(I,K)+W(I-1,K))*(V(I,K)+V(I,KM)-VB12)*RD4Z

c         y44(i)=-rho1(k)*am1(k)*r4dz2*y7(i)*(vv1(i,k)-vv1(i,km))
c         y66(i)=-rho1(k)*am1(k)*r4dz2*bsk4(k)*(vv1(i,k)-vv1(i,km)
c     1                                  -vb1(k)+vb1(km))

         y44(i)=rho1(k)*(-am1(k)*r4dz2*(y7(i)*(vv1(i,k)-vv1(i,km))
     1                   +bsk4(k)*(vv1(i,k)-vv1(i,km)-vb1(k)+vb1(km))))

c
        u1(i,km)=am(km)*rrho(km)*(y1(i)-y2(i)+y3(i)-y4(i))
        v1(i,km)=am(km)*rrho(km)*(y11(i)-y22(i)+y33(i)-y44(i))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         GRID_UZ(I,KM)=GRID_UZ(I,KM)+am(km)*rrho(km)*(Y1(I)-Y2(I))
c         DIFF_NU(I,KM)=DIFF_NU(I,KM)+am(km)*rrho(km)*(Y3(I)-Y4(I))
c         DIFF_UZ(I,KM)=DIFF_UZ(I,KM)+am(km)*rrho(km)*(Y5(I)-Y6(I))
C
c         GRID_VZ(I,KM)=GRID_VZ(I,KM)+am(km)*rrho(km)*(Y11(I)-Y22(I))
c         DIFF_NV(I,KM)=DIFF_NV(I,KM)+am(km)*rrho(km)*(Y33(I)-Y44(I))
c         DIFF_VZ(I,KM)=DIFF_VZ(I,KM)+am(km)*rrho(km)*(Y55(I)-Y66(I))
C
         GRID_UZ(I,KM)=am(km)*rrho(km)*(Y1(I)-Y2(I))
         DIFF_NU(I,KM)=am(km)*rrho(km)*(Y3(I)-Y4(I))
         DIFF_UZ(I,KM)=am(km)*rrho(km)*(Y5(I)-Y6(I))
C
         GRID_VZ(I,KM)=am(km)*rrho(km)*(Y11(I)-Y22(I))
         DIFF_NV(I,KM)=am(km)*rrho(km)*(Y33(I)-Y44(I))
         DIFF_VZ(I,KM)=am(km)*rrho(km)*(Y55(I)-Y66(I))
c
         y1(i)=y2(i)
         y3(i)=y4(i)
         y5(i)=y6(i)
         y11(i)=y22(i)
         y33(i)=y44(i)
         y55(i)=y66(i)
c
  200  continue
 2000 continue
c       a22k=am(kles)*rrho(kles)
      do 210 i=2,iles
       u1(i,kles)=am(kles)*rrho(kles)*(y2(i)+y4(i))
       v1(i,kles)=am(kles)*rrho(kles)*(y22(i)+y44(i))
C
c       GRID_UZ(I,KLES)=GRID_UZ(I,KLES)+am(kles)*rrho(kles)*Y2(i)
c       DIFF_NU(I,KLES)=DIFF_NU(I,KLES)+am(kles)*rrho(kles)*y4(i)
c       DIFF_UZ(I,KLES)=DIFF_UZ(I,KLES)+am(kles)*rrho(kles)*y6(i)
C
c       GRID_VZ(I,KLES)=GRID_VZ(I,KLES)+am(kles)*rrho(kles)*Y22(i)
c       DIFF_NV(I,KLES)=DIFF_NV(I,KLES)+am(kles)*rrho(kles)*y44(i)
c       DIFF_VZ(I,KLES)=DIFF_VZ(I,KLES)+am(kles)*rrho(kles)*y66(i)

       GRID_UZ(I,KLES)=am(kles)*rrho(kles)*Y2(i)
       DIFF_NU(I,KLES)=am(kles)*rrho(kles)*y4(i)
       DIFF_UZ(I,KLES)=am(kles)*rrho(kles)*y6(i)
C
       GRID_VZ(I,KLES)=am(kles)*rrho(kles)*Y22(i)
       DIFF_NV(I,KLES)=am(kles)*rrho(kles)*y44(i)
       DIFF_VZ(I,KLES)=am(kles)*rrho(kles)*y66(i)
C
  210 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  tao (11-13-97)
        if (ijkadv .eq. 0 .and. I_FOUR .EQ. 1) then
c
c
               r2dt2=0.005/dt
c
           do i=2,iles
             do k=2,kles
               y1(k)=uu1(i,k+1)-2.*uu1(i,k)+uu1(i,k-1)
     1               -ub1(k+1)+2.*ub1(k)-ub1(k-1)
               y2(k)=vv1(i,k+1)-2.*vv1(i,k)+vv1(i,k-1)
     1               -vb1(k+1)+2.*vb1(k)-vb1(k-1)
             enddo
             y1(1)=y1(2)
             y2(1)=y2(2)
             y1(kmax)=y1(kles)
             y2(kmax)=y2(kles)
             do k=2,kles
                u1(i,k)=u1(i,k)-r2dt2*(y1(k+1)-2.*y1(k)+y1(k-1))
                v1(i,k)=v1(i,k)-r2dt2*(y2(k+1)-2.*y2(k)+y2(k-1))
c
                DIFF_NU(I,K)=DIFF_NU(I,K)
     1                         -R2DT2*(Y1(K+1)-2.*Y1(K)+Y1(K-1))
                DIFF_NV(I,K)=DIFF_NV(I,K)
     1                         -R2DT2*(Y2(K+1)-2.*Y2(K)+Y2(K-1))
c
             enddo
           enddo
cc        IF (I_FOUR .EQ. 1) THEN
c
c               r2dt8=0.008/d2t
c               r2dt2=0.002/d2t
c
c           do i=2,iles
c             do k=3,kl2-1
c               y1(k)=(uu1(i,k+2)-3.*(uu1(i,k+1)-uu1(i,k))-uu1(i,k-1))
c     1               -(ub1(k+2)-3.*(ub1(k+1)-ub1(k))-ub1(k-1))
c               y2(k)=(vv1(i,k+2)-3.*(vv1(i,k+1)-vv1(i,k))-vv1(i,k-1))
c     1               -(vb1(k+2)-3.*(vb1(k+1)-vb1(k))-vb1(k-1))
c             enddo
c             do k=4,kl2-1
c                u1(i,k)=u1(i,k)+r2dt2*(y1(k-1)-y1(k))
c                v1(i,k)=v1(i,k)+r2dt2*(y2(k-1)-y2(k))
c             enddo
c             u1(i,3)=u1(i,3)+r2dt8*(uu1(i,4)-2.*uu1(i,3)+uu1(i,2))
c     1              -r2dt8*(ub1(4)-2.*ub1(3)+ub1(2))
c             v1(i,3)=v1(i,3)+r2dt8*(vv1(i,4)-2.*vv1(i,3)+vv1(i,2))
c             u1(i,kl2)=u1(i,kl2)+r2dt8*(uu1(i,kles)-2.*uu1(i,kl2)
c     1                          +uu1(i,kl2-1))
c     2                 -r2dt8*(ub1(kles)-2.*ub1(kl2)+ub1(kl2-1))
c
c             v1(i,kl2)=v1(i,kl2)+r2dt8*(vv1(i,kles)-2.*vv1(i,kl2)
c     1                          +vv1(i,kl2-1))
c     2                 -r2dt8*(vb1(kles)-2.*vb1(kl2)+vb1(kl2-1))
c          enddo

c
c        ELSE
c
c          y1(2)=0.
c          y1(kmax)=0.
c          y2(2)=0.
c          y2(kmax)=0.
c
c               r2dt8=0.008/d2t
c          do i=2,iles

c            do k=3,kles
c               y1(k)=uu1(i,k)-uu1(i,k-1)-ub1(k)+ub1(k-1)
c               y2(k)=vv1(i,k)-vv1(i,k-1)-vb1(k)+vb1(k-1)
c            enddo
c            do k=2,kles
c               u1(i,k)=u1(i,k)-r2dt8*(y1(k)-y1(k+1))
c               v1(i,k)=v1(i,k)-r2dt8*(y2(k)-y2(k+1))
c            enddo

c          enddo
C
C
c
cc        ENDIF
c
        endif
cccccc  tao (11-13-97)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nudge .ne. 1) then
        do 60 k=2,kles
          y1(k)=0.
   60     y3(k)=0.
        do 70 i=2,iles
        do 70 k=2,kles
          y1(k)=y1(k)+u1(i,k)
   70     y3(k)=y3(k)+v1(i,k)
        do 75 k=2,kles
          y1(k)=y1(k)*ril2
   75     y3(k)=y3(k)*ril2
        do 80 k=2,kles
        do 80 i=2,iles
          u1(i,k)=u1(i,k)-y1(k)
   80     v1(i,k)=v1(i,k)-y3(k)
      ENDIF
cc    ******   horizontal advections   *********************************
      do 2500 k=2,kles
       km=k-1
       kp=k+1
c
       y1(1)=a4k(k)*(uu1(3,k)-3.*(uu1(2,k)-uu1(1,k))-uu1(il2,k))
       y2(1)=a4k(k)*(vv1(3,k)-3.*(vv1(2,k)-vv1(1,k))-vv1(il2,k))
       do i=2,il2
         y1(i)=a4k(k)*(uu1(i+2,k)-3.*(uu1(i+1,k)-uu1(i,k))-uu1(i-1,k))
         y2(i)=a4k(k)*(vv1(i+2,k)-3.*(vv1(i+1,k)-vv1(i,k))-vv1(i-1,k))
       enddo
       do i=2,il2
         u1(i,k)=u1(i,k)+y1(i-1)-y1(i)
         v1(i,k)=v1(i,k)+y2(i-1)-y2(i)
c
          DIFF_NU(I,K)=DIFF_NU(I,K)+y1(i-1)-y1(i)
          DIFF_NV(I,K)=DIFF_NV(I,K)+y2(i-1)-y2(i)
c
       enddo
       u1(iles,k)=u1(iles,k)+y1(il2)-y1(1)
       v1(iles,k)=v1(iles,k)+y2(il2)-y2(1)
c
         DIFF_NU(ILES,K)=DIFF_NU(ILES,K)+Y1(IL2)-Y1(1)
         DIFF_NV(ILES,K)=DIFF_NV(ILES,K)+Y2(IL2)-Y2(1)
c
       do i=2,iles
          u1(i,k)=u1(i,k)-a_irf*rfa(k)*(uu1(i,k)-ub1(k))
     1                   +fcor*(v(i,k)-vb(k))
          v1(i,k)=v1(i,k)-a_irf*rfa(k)*(vv1(i,k)-vb1(k))
     2                   -fcor*(u(i,k)-ub(k))
c          u1(i,k)=u1(i,k)-a0000*am(k)*y0(i)
c     1          *(wb(kp)+wb(k))*(u(i,kp)-ub(kp)-u(i,km)+ub(km))*rd4z
c    2          -a0000*ub(k)*(u(i+1,k)-u(i-1,k))*rd2x
c          v1(i,k)=v1(i,k)-a0000*am(k)*y0(i)
c     1          *(wb(kp)+wb(k))*(v(i,kp)-vb(kp)-v(i,km)+vb(km))*rd4z
c     2          -a0000*ub(k)*(v(i+1,k)-v(i-1,k))*rd2x
c
          DIFF_NU(I,K)=DIFF_NU(I,K)-A_IRF*RFA(K)*(UU1(I,K)-UB1(K))
          DIFF_NV(I,K)=DIFF_NV(I,K)-A_IRF*RFA(K)*(VV1(I,K)-VB1(K))
c
       ENDDO
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       c33k=c3(k)*rdx
       ub2=0.0
c
       ub12=2.*ub(k)
       vb12=2.*vb(k)
c
       if (iadvh .eq. 2) then
       do 220 i=1,iles
c        y1(i)=(u(i,k)+u(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
c     1          -2.*ak(i,k)*(uu1(i+1,k)-uu1(i,k))*rdx2
c     2          +c33k*ak(i,k)*ak(i,k)
c        y2(i)=(v(i,k)+v(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
c     1          -(ak(i,k))*(vv1(i+1,k)-vv1(i,k))*rdx2
c
        y1(i)=(u(i,k)+u(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
        y3(i)=-2.*ak(i,k)*(uu1(i+1,k)-uu1(i,k))*rdx2
     1        +c33k*ak(i,k)*ak(i,k)
        y2(i)=(v(i,k)+v(i+1,k))*(u(i,k)+u(i+1,k)-ub2)*rd4x
        y4(i)=-(ak(i,k))*(vv1(i+1,k)-vv1(i,k))*rdx2
c
        y11(i)=(u(i,k)+u(i+1,k)-ub12)*(u(i,k)+u(i+1,k)-ub12)*rd4x
        y22(i)=(v(i,k)+v(i+1,k)-vb12)*(u(i,k)+u(i+1,k)-ub12)*rd4x
  220  continue
       do 230 i=2,iles
         u1(i,k)=u1(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i)
         v1(i,k)=v1(i,k)+y2(i-1)-y2(i)+y4(i-1)-y4(i)
c
c          GRID_UX(I,K)=GRID_UX(I,K)+Y1(I-1)-Y1(I)
c          DIFF_UX(I,K)=DIFF_UX(I,K)+Y11(I-1)-Y11(I)
c
c          GRID_VX(I,K)=GRID_VX(I,K)+Y2(I-1)-Y2(I)
c          DIFF_VX(I,K)=DIFF_VX(I,K)+Y22(I-1)-Y22(I)
c
c
          GRID_UX(I,K)=Y1(I-1)-Y1(I)
          DIFF_NU(I,K)=DIFF_NU(I,K)+Y3(I-1)-Y3(I)
          DIFF_UX(I,K)=Y11(I-1)-Y11(I)
c
          GRID_VX(I,K)=Y2(I-1)-Y2(I)
          DIFF_NV(I,K)=DIFF_NV(I,K)+Y4(I-1)-Y4(I)
          DIFF_VX(I,K)=Y22(I-1)-Y22(I)
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  230 continue
c
      else
c
cc    ******   4-th order horizontal advections   **********************
       do 240 i=1,iles
        y1(i)=-2.*ak(i,k)*(uu1(i+1,k)-uu1(i,k))*rdx2
     1         +c3(k)*ak(i,k)**2*rdx
        y2(i)=-ak(i,k)*(vv1(i+1,k)-vv1(i,k))*rdx2
  240  continue
       do 250 i=2,iles
        u1(i,k)=u1(i,k)+y1(i-1)-y1(i)
        v1(i,k)=v1(i,k)+y2(i-1)-y2(i)
c
          DIFF_NU(I,K)=DIFF_NU(I,K)+Y1(I-1)-Y1(I)
          DIFF_NV(I,K)=DIFF_NV(I,K)+Y2(I-1)-Y2(I)
c
  250  continue
       y1(1)=d516x*(u(1,k)+u(2,k))*(u(1,k)+u(2,k)-ub2)
     1       -d32x*(u(il2,k)+u(3,k))*(u(il2,k)+u(3,k)-ub2)
       y2(1)=d516x*(v(1,k)+v(2,k))*(u(1,k)+u(2,k)-ub2)
     1       -d32x*(v(il2,k)+v(3,k))*(u(il2,k)+u(3,k)-ub2)
c
       y11(1)=d516x*(u(1,k)+u(2,k)-ub12)*(u(1,k)+u(2,k)-ub12)
     1       -d32x*(u(il2,k)+u(3,k)-ub12)*(u(il2,k)+u(3,k)-ub12)
       y22(1)=d516x*(v(1,k)+v(2,k)-vb12)*(u(1,k)+u(2,k)-ub12)
     1       -d32x*(v(il2,k)+v(3,k)-vb12)*(u(il2,k)+u(3,k)-ub12)
c
       do 260 i=2,il2
        y1(i)=d516x*(u(i,k)+u(i+1,k))*(u(i,k)+u(i+1,k)-ub2)
     1        -d32x*(u(i-1,k)+u(i+2,k))*(u(i-1,k)+u(i+2,k)-ub2)
        y2(i)=d516x*(v(i,k)+v(i+1,k))*(u(i,k)+u(i+1,k)-ub2)
     1        -d32x*(v(i-1,k)+v(i+2,k))*(u(i-1,k)+u(i+2,k)-ub2)
c
        y11(i)=d516x*(u(i,k)+u(i+1,k)-ub12)*(u(i,k)+u(i+1,k)-ub12)
     1        -d32x*(u(i-1,k)+u(i+2,k)-ub12)*(u(i-1,k)+u(i+2,k)-ub12)
        y22(i)=d516x*(v(i,k)+v(i+1,k)-vb12)*(u(i,k)+u(i+1,k)-ub12)
     1        -d32x*(v(i-1,k)+v(i+2,k)-vb12)*(u(i-1,k)+u(i+2,k)-ub12)
  260  continue
       y3(1)=-d96x*(u(il2,k)+u(1,k))*(u(il2,k)+u(1,k)-ub2)
       y4(1)=-d96x*(v(il2,k)+v(1,k))*(u(il2,k)+u(1,k)-ub2)
c
       y33(1)=-d96x*(u(il2,k)+u(1,k)-ub12)*(u(il2,k)+u(1,k)-ub12)
       y44(1)=-d96x*(v(il2,k)+v(1,k)-vb12)*(u(il2,k)+u(1,k)-ub12)
       do 270 i=2,imax
        y3(i)=-d96x*(u(i-1,k)+u(i,k))*(u(i-1,k)+u(i,k)-ub2)
        y4(i)=-d96x*(v(i-1,k)+v(i,k))*(u(i-1,k)+u(i,k)-ub2)
c
        y33(i)=-d96x*(u(i-1,k)+u(i,k)-ub12)*(u(i-1,k)+u(i,k)-ub12)
        y44(i)=-d96x*(v(i-1,k)+v(i,k)-vb12)*(u(i-1,k)+u(i,k)-ub12)
  270  continue
       do 280 i=2,il2
        u1(i,k)=u1(i,k)+y1(i-1)-y1(i)+y3(i-1)-y3(i+2)
        v1(i,k)=v1(i,k)+y2(i-1)-y2(i)+y4(i-1)-y4(i+2)
c
c          GRID_UX(I,K)=GRID_UX(I,K)+Y1(I-1)-Y1(I)+Y3(I-1)-Y3(I+2)
c          GRID_VX(I,K)=GRID_VX(I,K)+Y2(I-1)-Y2(I)+Y4(I-1)-Y4(I+2)
c
c          DIFF_UX(I,K)=DIFF_UX(I,K)+Y11(I-1)-Y11(I)+Y33(I-1)-Y33(I+2)
c          DIFF_VX(I,K)=DIFF_VX(I,K)+Y22(I-1)-Y22(I)+Y44(I-1)-Y44(I+2)
c
          GRID_UX(I,K)=Y1(I-1)-Y1(I)+Y3(I-1)-Y3(I+2)
          GRID_VX(I,K)=Y2(I-1)-Y2(I)+Y4(I-1)-Y4(I+2)
c
          DIFF_UX(I,K)=Y11(I-1)-Y11(I)+Y33(I-1)-Y33(I+2)
          DIFF_VX(I,K)=Y22(I-1)-Y22(I)+Y44(I-1)-Y44(I+2)
  280  continue
       u1(iles,k)=u1(iles,k)+y1(il2)-y1(1)+y3(il2)-y3(3)
       v1(iles,k)=v1(iles,k)+y2(il2)-y2(1)+y4(il2)-y4(3)
c
c        GRID_UX(ILES,K)=GRID_UX(ILES,K)+Y1(IL2)-Y1(1)+Y3(IL2)-Y3(3)
c        GRID_VX(ILES,K)=GRID_VX(ILES,K)+Y2(IL2)-Y2(1)+Y4(IL2)-Y4(3)
c
c        DIFF_UX(ILES,K)=DIFF_UX(ILES,K)+Y11(IL2)-Y11(1)+Y33(IL2)-Y33(3)
c        DIFF_VX(ILES,K)=DIFF_VX(ILES,K)+Y22(IL2)-Y22(1)+Y44(IL2)-Y44(3)
c
        GRID_UX(ILES,K)=Y1(IL2)-Y1(1)+Y3(IL2)-Y3(3)
        GRID_VX(ILES,K)=Y2(IL2)-Y2(1)+Y4(IL2)-Y4(3)
c
        DIFF_UX(ILES,K)=Y11(IL2)-Y11(1)+Y33(IL2)-Y33(3)
        DIFF_VX(ILES,K)=Y22(IL2)-Y22(1)+Y44(IL2)-Y44(3)
       endif
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 2500 continue
c
c
ctao 9/22/97
       if (nudge .ne. 1) then
         do 295 k=2,kles
          do i=2,iles
c
c             U_LARGE(I,K)=U_LARGE(I,K)+UBT(K)
c             V_LARGE(I,K)=V_LARGE(I,K)+VBT(K)
             U_LARGE(I,K)=UBT(K)
             V_LARGE(I,K)=VBT(K)
c
             u1(i,k)=u1(i,k)+ubt(k)
             v1(i,k)=v1(i,k)+vbt(k)
          enddo
  295    continue
       endif
        if (nudge .eq. 1) then
          do k=2,kles
          do i=2,iles
c
c             U_LARGE(I,K)=U_LARGE(I,K)-FACTOR_NUDING*(UB(K)-UBI(K))
c             V_LARGE(I,K)=V_LARGE(I,K)-FACTOR_NUDING*(VB(K)-VBI(K))
c
            u1(i,k)=u1(i,k)-factor_nuding*(ub(k)-ubi(k))
            v1(i,k)=v1(i,k)-factor_nuding*(vb(k)-vbi(k))
c
             U_LARGE(I,K)=-FACTOR_NUDING*(UB(K)-UBI(K))
             V_LARGE(I,K)=-FACTOR_NUDING*(VB(K)-VBI(K))
         enddo
         enddo
       endif
ctao 9/22/97
c
c
cc    ******************************************************************
      do 300 k=2,kles
      do 300 i=2,iles
         VTEM_DT(i,k)=vv1(i,k)
       y1(i)=v(i,k)+eps*(-2.*v(i,k)+vv1(i,k))
       v(i,k)=vv1(i,k)+v1(i,k)*d2t
       vv1(i,k)=y1(i)+eps*v(i,k)
  300 continue
      do 400 k=1,kmax
      do 400 i=1,imax
  400   ak(i,k)=2.*ak(i,k)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine advectak (x,x1)

      implicit none
CC    ****   COMPUTE ADVECTION OF K COEFFICIENT
      integer nx,nz
      PARAMETER (NX=514,NZ=43)

      real    x(nx,nz),x1(nx,nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    u1(nx,nz),w1(nx,nz)
      common/bsat/ u1
      common/bsat1/ w1

      real    uu1(nx,nz),ww1(nx,nz)
      common/b2u/ uu1
      common/b2w/ ww1

      real    umd(nx,nz),vmd(nx,nz),wmd(nx,nz),umod(nx,nz),wmod(nx,nz)
      common/dumuw/ umd,vmd,wmd
      common/bw0/ umod,wmod

      real    dxxt(nx),dzzt(nz),dxr(nx),dzr(nz)
      common/bcor2/ dxxt,dzzt,dxr,dzr

      real    tb(nz),qb(nz),rho1(nz),rho(nz),zz(nz),zz1(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,zz,zz1,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k

      save


CC    **************************************************************
        DO 10 I=1,IMAX
         DXR(I)=RDX
   10    DXXT(I)=DT*RDX
       DO 20 K=2,KLES
        DZR(K)=AM(K)*RDZ/RHO(K)
   20   DZZT(K)=DT*AM1(K)*RDZ
      DO 100 K=1,KMAX
      DO 100 I=1,IMAX
        UMOD(I,K)=0.
        WMOD(I,K)=0.
       U1(I,K)=.5*(3.*UU1(I,K)-UMD(I,K))
       W1(I,K)=.5*(3.*WW1(I,K)-WMD(I,K))
  100  X1(I,K)=X(I,K)
      DO 150 K=1,KMAX
      DO 150 I=1,IMAX
        IF (X(I,K) .LT. 1.E-5) X(I,K)=0.0
       X1(I,K)=X(I,K)
  150 CONTINUE
      CALL FADVAK (X,U1,W1)
       CALL BOUNDY (X,X1)
       CALL FADVUW (X,X1,U1,W1,UMOD,WMOD)
       CALL FADVAK (X,UMOD,WMOD)
       CALL BOUNDY (X,X1)
      DO 200 K=1,KMAX
      DO 200 I=1,IMAX
        IF (X(I,K) .LT. 1.E-5) X(I,K)=0.
  200 CONTINUE
      RETURN
      END
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmpf
      implicit none
c     ******   compute force terms   ***********************************
c     ******   update w, u and v without pressure force

      integer nx,nz,nx10
      parameter (NX=514,NZ=43)
      parameter (nx10=10*nx)

      integer LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    uu1(nx,nz),vv1(nx,nz),ww1(nx,nz),
     $        u  (nx,nz),v  (nx,nz),w  (nx,nz)
      common/b1u/ u
      common/b1v/ v
      common/b1w/ w
      common/b2u/ uu1
      common/b2v/ vv1
      common/b2w/ ww1

      real    f(nx,nz),w1(nx,nz)
      common/bsat/ f
      common/bsat1/ w1

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx10)
      common/ba/ y1,y2,y3,y4,y5,y6,y7

      real UTEM_DT(NX,NZ),VTEM_DT(NX,NZ)
      COMMON/DEBUG_UUVV/ UTEM_DT,VTEM_DT

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a1,a2,rd2t

      save

c     ********************************
      rd2t=1./d2t
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO 10 K=2,KLES
        Y1(K)=AM(K)*RDZ*RRHO(K)
        IF (LIPPS .eq. 0) THEN
          Y2(K)=DZ2/(CP*TB(K))
        else
          y2(k)=dz2/cp
        endif
   10   F(IMAX,K)=F(2,K)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 1000 k=2,kles
        do 20 i=2,imax
   20    y3(i)=f(i,k)
       do 100 i=2,iles
       f(i,k)=((y3(i+1)-y3(i))*rdx
     1        +(w1(i,k+1)*rho1(k+1)-w1(i,k)*rho1(k))*y1(k)
     2       +((uu1(i+1,k)-uu1(i,k))*rdx
     3       +(ww1(i,k+1)*rho1(k+1)-ww1(i,k)*rho1(k))*y1(k))*rd2t)*y2(k)

          UTEM_DT(I,K)=UU1(I,K)

         y4(i)=u(i,k)+eps*(-2.*u(i,k)+uu1(i,k))
        u(i,k)=uu1(i,k)+y3(i)*d2t
        uu1(i,k)=y4(i)
c
         y5(i)=w(i,k)+eps*(-2.*w(i,k)+ww1(i,k))
        w(i,k)=ww1(i,k)+w1(i,k)*d2t
  100   ww1(i,k)=y5(i)
 1000 continue
c   ****   set boundary conditions  ************************************
      IF (LIPPS .eq. 0) THEN
        a1=dz*ba(2)/(cp*tb(2)*am1(2))
        a2=dz*bb(kles)/(cp*tb(kles)*am1(kles))
      ELSE
        A1=DZ*BA(2)/(CP*AM1(2))
        A2=DZ*BB(KLES)/(CP*AM1(KLES))
      ENDIF
       do 60 i=2,iles
        f(i,2)=f(i,2)+w1(i,2)*a1
   60   f(i,kles)=f(i,kles)-w1(i,kmax)*a2
C
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cmpwuv (id)

      implicit none
c     *******   compute new w, u and v   *****************************

      integer nx,nz,nx10,itt,nt,id,mt
c     parameter (NX=514,NZ=43,ITT=244)
      parameter (NX=514,NZ=43,NT=2880,ITT=244) !6/21/01
      parameter (nx10=10*nx)

      integer LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    uu1(nx,nz),vv1(nx,nz),ww1(nx,nz),
     $        u  (nx,nz),v  (nx,nz),w  (nx,nz)
      common/b1u/ u
      common/b1v/ v
      common/b1w/ w
      common/b2u/ uu1
      common/b2v/ vv1
      common/b2w/ ww1

      real    qcl(nx,nz),qrn(nx,nz),qci(nx,nz),qcs(nx,nz),qcg(nx,nz)
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1i/ qci
      common/b1s/ qcs
      common/b1g/ qcg

      real    pi(nx,nz)
      common/bsat/ pi

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nx10)
      common/ba/ y1,y2,y3,y4,y5,y6,y7

      real ppress(nx,nz)
      COMMON/SUE/ PPRESS

      real    tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/bb6/ tls1,tls2,qls1,qls2,tls3,tls4,qls3,qls4,sft,sfq,wbt,
     $            wb_6h,ub_6h,ubt,q1_6h,q1t,q2_6h,q2t,vb_6h,vbt
c
      integer     IT(NX),IV(NT),ICS(NX,4),IBZ(NX,4)
      COMMON/BCH/ IT,IV,ICS,IBZ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ICS(I,3)=1:  STRATIFORM REGION                                  C
C      ICS(I,2)=1:  CONVECTIVE REGION                                  C
C      ICS(I,4)=1:  NO SFC RAIN BUT CLOUD ALOFT                        C
C      ICE(I,1)=1:  TOTAL MODEL DOMAIN                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real            DIFF_VX(NX,NZ),DIFF_VZ(NX,NZ),GRID_VX(NX,NZ),
     1                GRID_VZ(NX,NZ),DIFF_NV(NX,NZ),V_LARGE(NX,NZ),
     2                PRE_V(NX,NZ),DT_VWIND(NX,NZ)
      COMMON/DEBUG_V/ DIFF_VX,DIFF_VZ,GRID_VX,
     1                GRID_VZ,DIFF_NV,V_LARGE,
     2                PRE_V,DT_VWIND

      real            DIFF_UX(NX,NZ),DIFF_UZ(NX,NZ),GRID_UX(NX,NZ),
     1                GRID_UZ(NX,NZ),DIFF_NU(NX,NZ),U_LARGE(NX,NZ),
     2                PRE_U(NX,NZ),DT_UWIND(NX,NZ)
      COMMON/DEBUG_U/ DIFF_UX,DIFF_UZ,GRID_UX,
     1                GRID_UZ,DIFF_NU,U_LARGE,
     2                PRE_U,DT_UWIND

      real             SDIFF_VX(NZ,4),SDIFF_VZ(NZ,4),SGRID_VX(NZ,4),
     1                 SGRID_VZ(NZ,4),SDIFF_NV(NZ,4),SV_LARGE(NZ,4),
     2                 SPRE_V(NZ,4),SDT_V(NZ,4),TN_VWIND(NZ,4)
      COMMON/DEBUG_SV/ SDIFF_VX,SDIFF_VZ,SGRID_VX,
     1                 SGRID_VZ,SDIFF_NV,SV_LARGE,
     2                 SPRE_V,SDT_V,TN_VWIND

      real             SDIFF_UX(NZ,4),SDIFF_UZ(NZ,4),SGRID_UX(NZ,4),
     1                 SGRID_UZ(NZ,4),SDIFF_NU(NZ,4),SU_LARGE(NZ,4),
     2                 SPRE_U(NZ,4),SDT_U(NZ,4),TN_UWIND(NZ,4)
      COMMON/DEBUG_SU/ SDIFF_UX,SDIFF_UZ,SGRID_UX,
     1                 SGRID_UZ,SDIFF_NU,SU_LARGE,
     2                 SPRE_U,SDT_U,TN_UWIND

      real        UB_MEAN(NZ,4),VB_MEAN(NZ,4)

      COMMON/SUV/ UB_MEAN,VB_MEAN

      real               UTEM_DT(NX,NZ),VTEM_DT(NX,NZ)
      COMMON/DEBUG_UUVV/ UTEM_DT,VTEM_DT

      real             SSDIFF_VX(NZ,4),SSDIFF_VZ(NZ,4),SSGRID_VX(NZ,4),
     1                 SSGRID_VZ(NZ,4),SSDIFF_NV(NZ,4),SSV_LARGE(NZ,4),
     2                 SSPRE_V(NZ,4),SSDT_V(NZ,4),STN_VWIND(NZ,4)
      COMMON/DEBUG_SVS/ SSDIFF_VX,SSDIFF_VZ,SSGRID_VX,
     1                  SSGRID_VZ,SSDIFF_NV,SSV_LARGE,
     2                  SSPRE_V,SSDT_V,STN_VWIND

      real             SSDIFF_UX(NZ,4),SSDIFF_UZ(NZ,4),SSGRID_UX(NZ,4),
     1                 SSGRID_UZ(NZ,4),SSDIFF_NU(NZ,4),SSU_LARGE(NZ,4),
     2                 SSPRE_U(NZ,4),SSDT_U(NZ,4),STN_UWIND(NZ,4)
      COMMON/DEBUG_SUS/ SSDIFF_UX,SSDIFF_UZ,SSGRID_UX,
     1                  SSGRID_UZ,SSDIFF_NU,SSU_LARGE,
     2                  SSPRE_U,SSDT_U,STN_UWIND

      real         SUB_MEAN(NZ,4),SVB_MEAN(NZ,4)

      COMMON/SSUV/ SUB_MEAN,SVB_MEAN

      real qc_t,qc_tl,ww_t

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    a1,a2,a11,a22,a33,a44,a55,a66,a77,a88,a99,ab11
      real    b11,b22,b33,b44,b55,b66,b77,b88,b99,ab22
      save
c     ******************************************************************
        a1=.5*cp*d2t*rdz
        a2=cp*d2t*rdx
        IF (LIPPS .EQ. 1) A1=CP*D2T*RDZ
       do 10 k=2,kles
        y1(k)=am1(k)*(tb(k-1)+tb(k))*a1
        y2(k)=tb(k)*a2
        IF (LIPPS .EQ. 1) then
          Y1(K)=AM1(K)*A1
          y2(k)=a2
        endif
   10   pi(1,k)=pi(iles,k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 1000 k=2,kles
       if(k.eq.2) go to 150
       do 100 i=2,iles
        w(i,k)=w(i,k)-y1(k)*(pi(i,k)-pi(i,k-1))
  100   ww1(i,k)=ww1(i,k)+eps*w(i,k)
  150  continue
       do 175 i=2,iles
          PPRESS(I,K)=PI(I,K)
c
c         PRE_U(I,K)=PRE_U(I,K)-Y2(K)*(PI(I,K)-PI(I-1,K))/d2t
         PRE_U(I,K)=-Y2(K)*(PI(I,K)-PI(I-1,K))/d2t
c
        u(i,k)=u(i,k)-y2(k)*(pi(i,k)-pi(i-1,k))
  175   uu1(i,k)=uu1(i,k)+eps*u(i,k)
 1000 continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      if (iuvbar .eq. 1) then
c        do k=2,kles
c        do i=2,iles
c          u(i,k)=u(i,k)+ubt(k)*dt
c          uu1(i,k)=uu1(i,k)+ubt(k)*dt
c          v(i,k)=v(i,k)+vbt(k)*dt
c          vv1(i,k)=vv1(i,k)+vbt(k)*dt
c        enddo
c        enddo
c      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 20 i=2,iles
       w(i,2)=0.
       ww1(i,2)=0.
   20 continue
      call boundy (u,uu1)
      call boundy (v,vv1)
      call boundy (w,ww1)
      do 30 i=1,imax
       w(i,2)=0.
       ww1(i,2)=0.
       w(i,kmax)=0.
       ww1(i,kmax)=0.
   30 continue
cc    ***********************
      do 40 k=1,kmax
      y1(k)=0.
      y2(k)=0.
      y3(k)=0.
   40 y4(k)=0.
      do 50 i=2,iles
      do 50 k=1,kmax
       y1(k)=y1(k)+u(i,k)
       y2(k)=y2(k)+uu1(i,k)
       y3(k)=y3(k)+v(i,k)
   50  y4(k)=y4(k)+vv1(i,k)
      do 60 k=1,kmax
       ub(k)=y1(k)*ril2
       ub1(k)=y2(k)*ril2
       vb(k)=y3(k)*ril2
       vb1(k)=y4(k)*ril2
   60 continue
c
      DO K=2,KLES
         DO I=2,ILES
c           DT_UWIND(I,K)=DT_UWIND(I,K)+(U(I,K)-UU1(I,K))/DT
c           DT_VWIND(I,K)=DT_VWIND(I,K)+(V(I,K)-VV1(I,K))/DT
           DT_UWIND(I,K)=(UU1(I,K)-UTEM_DT(I,K))/D2T
           DT_VWIND(I,K)=(VV1(I,K)-VTEM_DT(I,K))/D2T
         ENDDO
      ENDDO
C
CC
C
      IF (ID .EQ. 1) THEN

          do i=2,iles
             ICS(I,1)=1
          enddo

         DO MT=1,4

            DO K=2,KLES
            A11=0.
            A22=0.
            A33=0.
            A44=0.
            A55=0.
            A66=0.
            A77=0.
            A88=0.
            A99=0.
c
            B11=0.
            B22=0.
            B33=0.
            B44=0.
            B55=0.
            B66=0.
            B77=0.
            B88=0.
            B99=0.
c
            AB11=0.
            AB22=0.
               DO I=2,ILES
                  IF (ICS(I,MT) .EQ. 1) THEN
                     A11=A11+DIFF_VX(I,K)
                     A22=A22+DIFF_VZ(I,K)
                     A33=A33+GRID_VX(I,K)
                     A44=A44+GRID_VZ(I,K)
                     A55=A55+DIFF_NV(I,K)
                     A66=A66+V_LARGE(I,K)
                     A77=A77+PRE_V(I,K)
                     A88=A88+DT_VWIND(I,K)
                     A99=A99+1.
                     AB11=AB11+V(I,K)
CCCCC
                     B11=B11+DIFF_UX(I,K)
                     B22=B22+DIFF_UZ(I,K)
                     B33=B33+GRID_UX(I,K)
                     B44=B44+GRID_UZ(I,K)
                     B55=B55+DIFF_NU(I,K)
                     B66=B66+U_LARGE(I,K)
                     B77=B77+PRE_U(I,K)
                     B88=B88+DT_UWIND(I,K)
                     B99=B99+1.
                     AB22=AB22+U(I,K)
                  ENDIF
               ENDDO
C
               SDIFF_VX(K,MT)=SDIFF_VX(K,MT)+A11
               SDIFF_VZ(K,MT)=SDIFF_VZ(K,MT)+A22
               SGRID_VX(K,MT)=SGRID_VX(K,MT)+A33
               SGRID_VZ(K,MT)=SGRID_VZ(K,MT)+A44
               SDIFF_NV(K,MT)=SDIFF_NV(K,MT)+A55
               SV_LARGE(K,MT)=SV_LARGE(K,MT)+A66
               SPRE_V(K,MT)=SPRE_V(K,MT)+A77
               SDT_V(K,MT)=SDT_V(K,MT)+A88
               TN_VWIND(K,MT)=TN_VWIND(K,MT)+A99
c 7/27/01 tao, collect VB_MEAN, UB_MEAN every "id", 120 sec= 2 min
c                IF (A99 .EQ. 0.) A99=1.E20 ! tao 7/27/01
c              VB_MEAN(K,MT)=VB_MEAN(K,MT)+AB11/A99 ! tao 7/27/01, every "id",  120 sec= 2 min
c 8/4/01 tao
           IF (A99 .GE. 0.01) VB_MEAN(K,MT)=VB_MEAN(K,MT)+AB11/A99
CCCCC
               SDIFF_UX(K,MT)=SDIFF_UX(K,MT)+B11
               SDIFF_UZ(K,MT)=SDIFF_UZ(K,MT)+B22
               SGRID_UX(K,MT)=SGRID_UX(K,MT)+B33
               SGRID_UZ(K,MT)=SGRID_UZ(K,MT)+B44
               SDIFF_NU(K,MT)=SDIFF_NU(K,MT)+B55
               SU_LARGE(K,MT)=SU_LARGE(K,MT)+B66
               SPRE_U(K,MT)=SPRE_U(K,MT)+B77
               SDT_U(K,MT)=SDT_U(K,MT)+B88
               TN_UWIND(K,MT)=TN_UWIND(K,MT)+B99
           IF (B99 .GE. 0.01) UB_MEAN(K,MT)=UB_MEAN(K,MT)+AB22/B99
            ENDDO
         ENDDO
C
cc
c
         do mt=1,4
            do i=1,nx
c              ibx(i,mt)=0
               ibz(i,mt)=0 ! by shie 8/1/01
            enddo
         enddo
          DO K=2,KLES
            do i=2,iles
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            qc_t=qcl(i,k)+qrn(i,k)+qci(i,k)+qcs(i,k)+qcg(i,k)
     1          +qcl(i-1,k)+qrn(i-1,k)+qci(i-1,k)+qcs(i-1,k)+qcg(i-1,k) ! shie
c    1            +qcl(i-1,k)+qrn(i-1,k)+qci(i-1,k)+qcs(i-1,k)+qcg(i-1,k)
c 8/6/01 tao, shie qc_tl (saturated) and qc_t (all) were found no significant
c         difference in impact, so comment out the sorting for qc_tl.
c           qc_tl=qcl(i,k)+qci(i,k)+qcl(i-1,k)+qci(i-1,k)
            ww_t=.0025*(w(i,k)+w(i,k-1)+w(i-1,k)+w(i-1,k-1))

              if (qc_t .ge. 2.e-5 .and. ww_t .ge. 1.) ibz(i,1)=1
              if (qc_t .ge. 2.e-5 .and. ww_t .le. -0.5) ibz(i,2)=1
              if (qc_t .ge. 2.e-5 .and. ww_t .ge. .001) ibz(i,3)=1
              if (qc_t .ge. 2.e-5 .and. ww_t .le. -0.001) ibz(i,4)=1
c             if (qc_tl .ge. 2.e-5 .and. ww_t .ge. 1.) ibz(i,3)=1
c             if (qc_tl .ge. 2.e-5 .and. ww_t .le. -0.5) ibz(i,4)=1
            enddo
           ENDDO
c
              DO MT=1,4
c
         DO K=2,KLES
c
            A11=0.
            A22=0.
            A33=0.
            A44=0.
            A55=0.
            A66=0.
            A77=0.
            A88=0.
            A99=0.
c
            B11=0.
            B22=0.
            B33=0.
            B44=0.
            B55=0.
            B66=0.
            B77=0.
            B88=0.
            B99=0.
c
            AB11=0.
            AB22=0.
              do i=2,iles
                IF (IBZ(I,MT) .EQ. 1) THEN
                   A11=A11+DIFF_VX(I,K)
                   A22=A22+DIFF_VZ(I,K)
                   A33=A33+GRID_VX(I,K)
                   A44=A44+GRID_VZ(I,K)
                   A55=A55+DIFF_NV(I,K)
                   A66=A66+V_LARGE(I,K)
                   A77=A77+PRE_V(I,K)
                   A88=A88+DT_VWIND(I,K)
                   A99=A99+1.
                   AB11=AB11+V(I,K)
CCCCC
                   B11=B11+DIFF_UX(I,K)
                   B22=B22+DIFF_UZ(I,K)
                   B33=B33+GRID_UX(I,K)
                   B44=B44+GRID_UZ(I,K)
                   B55=B55+DIFF_NU(I,K)
                   B66=B66+U_LARGE(I,K)
                   B77=B77+PRE_U(I,K)
                   B88=B88+DT_UWIND(I,K)
                   B99=B99+1.
                   AB22=AB22+U(I,K)
                ENDIF
              ENDDO
C
               SSDIFF_VX(K,MT)=SSDIFF_VX(K,MT)+A11
               SSDIFF_VZ(K,MT)=SSDIFF_VZ(K,MT)+A22
               SSGRID_VX(K,MT)=SSGRID_VX(K,MT)+A33
               SSGRID_VZ(K,MT)=SSGRID_VZ(K,MT)+A44
               SSDIFF_NV(K,MT)=SSDIFF_NV(K,MT)+A55
               SSV_LARGE(K,MT)=SSV_LARGE(K,MT)+A66
               SSPRE_V(K,MT)=SSPRE_V(K,MT)+A77
               SSDT_V(K,MT)=SSDT_V(K,MT)+A88
               STN_VWIND(K,MT)=STN_VWIND(K,MT)+A99
c                IF (A99 .EQ. 0.) A99=1.E20
c              SVB_MEAN(K,MT)=SVB_MEAN(K,MT)+AB11/A99
         IF (A99 .GE. 0.01) SVB_MEAN(K,MT)=SVB_MEAN(K,MT)+AB11/A99 ! shie
CCCCC
               SSDIFF_UX(K,MT)=SSDIFF_UX(K,MT)+B11
               SSDIFF_UZ(K,MT)=SSDIFF_UZ(K,MT)+B22
               SSGRID_UX(K,MT)=SSGRID_UX(K,MT)+B33
               SSGRID_UZ(K,MT)=SSGRID_UZ(K,MT)+B44
               SSDIFF_NU(K,MT)=SSDIFF_NU(K,MT)+B55
               SSU_LARGE(K,MT)=SSU_LARGE(K,MT)+B66
               SSPRE_U(K,MT)=SSPRE_U(K,MT)+B77
               SSDT_U(K,MT)=SSDT_U(K,MT)+B88
               STN_UWIND(K,MT)=STN_UWIND(K,MT)+B99
c                IF (B99 .EQ. 0.) B99=1.E20
c              SUB_MEAN(K,MT)=SUB_MEAN(K,MT)+AB22/B99
         IF (B99 .GE. 0.01) SUB_MEAN(K,MT)=SUB_MEAN(K,MT)+AB22/B99 ! shie
           ENDDO
           ENDDO

      ENDIF

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine domain (x,i1,i2,k1,k2)

      implicit none

      integer nx,nz
      parameter (NX=514,NZ=43)

      integer i1,i2,k1,k2
      real    x(nx,nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    y

      save

      i1=0
      y=0.
      do 100 i=2,iles
       do 10 k=2,kles
   10   y=max(y,x(i,k))
       if(y .gt. 0.) then
        i1=i
        go to 15
       endif
  100 continue
   15 if(i1 .eq. 0) then
       i1=2
       i2=2
       k1=2
       k2=2
       return
      endif
      y=0.
      do 200 i=iles,2,-1
       do 20 k=2,kles
   20  y=max(y,x(i,k))
       if(y .gt. 0.) then
        i2=i
        go to 25
       endif
  200 continue
   25 y=0.
      do 300 k=2,kles
       do 30 i=2,iles
   30  y=max(y,x(i,k))
       if(y .gt. 0.) then
        k1=k
        go to 35
       endif
  300 continue
   35 y=0.
      do 400 k=kles,2,-1
       do 40 i=2,iles
   40  y=max(y,x(i,k))
       if(y .gt. 0.) then
        k2=k
        go to 45
       endif
  400 continue
   45 i1=max(2,i1-3)
      i2=min(iles,i2+3)
      k1=max(2,k1-3)
      k2=min(kles,k2+3)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine slvpi (iflg)

      implicit none
c     ******   solve 2-d pressure equation  (fast fourier transform)

      integer nx,nz,nxm,nzm,nx12,iflg
      parameter (NX=514,NZ=43)
      parameter (nxm=nx-1,nzm=nz-1,nx12=12*nx)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    dt,d2t,ril2,f5(11)
      common/bb/ dt,d2t,ril2,f5

      real    f(nx,nz),ee(nx,nz)
      common/bsat/ f
      common/bsat1/ ee

      real    aux(nx,nz),ff(nx,nz)
      common/badv/ aux
      common/badv1/ ff

      real    tb(nz),qb(nz),rho1(nz),rho(nz),cc(nz),aa(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,cc,aa,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx12)
      common/ba/ y1,y2,y3,y4,y5

      real    ar(nxm,nzm),ai(nxm,nzm)
      common/rfft/ ar,ai

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    beta(nx)
      integer i,ihf,ihf1,ihf2,imi,k,k1
      real    pix

      save

c     ******   initialize the para for solver
      if(iflg.eq.1) go to 100
      ihf=il2/2+1
      ihf1=ihf+1
      ihf2=ihf1+1
c     ******
      pix=8.*atan(1.)*ril2
       beta(1)=0.
      do 20 i=2,iles
   20  beta(i)=2.*dz2*(cos(pix*(i-2))-1.)*rdx2
       beta(imax)=0.
      call fft (0,il2)
      return
c     ****   forward transformation in x
  100 continue
      do 22 k=2,kles
      do 22 i=2,iles
       ar(i-1,k-1)=f(i,k)*ril2
       ai(i-1,k-1)=0.
   22 continue
      call fft (1,il2)
      do 200 k=2,kles
       do 24 i=2,ihf1
   24   f(i,k)=ar(i-1,k-1)
       do 26 i=ihf2,iles
   26   f(i,k)=ai(i-1,k-1)
  200 continue
c     ******      ******
      do 30 i=2,iles
       aux(i,2)=1./(aa(2)-beta(i))
       ee(i,2)=aa(2)*aux(i,2)
       ff(i,2)=-f(i,2)*aux(i,2)
   30 continue
      do 300 k=3,kl2
      do 300 i=2,iles
       aux(i,k)=1./(aa(k)+cc(k)-beta(i)-cc(k)*ee(i,k-1))
       ee(i,k)=aa(k)*aux(i,k)
  300  ff(i,k)=(-f(i,k)+cc(k)*ff(i,k-1))*aux(i,k)
      y1(2)=0.
      do 32 i=3,iles
       y1(i)=(-f(i,kles)+cc(kles)*ff(i,kl2))
     1       /((1.-ee(i,kl2))*cc(kles)-beta(i))
   32 continue
      do 34 i=2,iles
   34  f(i,kles)=y1(i)
      do 350 k=3,kles
       k1=kmax+1-k
      do 350 i=2,iles
       y1(i)=ee(i,k1)*y1(i)+ff(i,k1)
  350  f(i,k1)=y1(i)
c     ****   backward transformation in x    ******
      do 400 k=2,kles
       do 40 i=2,ihf1
   40   ar(i-1,k-1)=f(i,k)
       do 42 i=ihf2,iles
   42   ai(i-1,k-1)=-f(i,k)
  400 continue
      do 44 k=1,kl2
       ai(1,k)=0.
       ai(ihf,k)=0.
   44 continue
      do 440 i=ihf1,il2
       imi=imax-i
      do 440 k=1,kl2
       ar(i,k)=ar(imi,k)
  440  ai(imi,k)=-ai(i,k)
      call fft (1,il2)
      do 500 k=2,kles
      do 500 i=2,iles
  500  f(i,k)=ar(i-1,k-1)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fft (iflg,l)

      implicit none

      integer nx,nz,nxi,nxm,nzm,nxq,nb,iflg,l
      parameter (NX=514,NZ=43)
      parameter (nxi=nx-2,nxm=nx-1,nzm=nz-1,nxq=nxi/2+1,nb=10*nx)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),y6(nx),y7(nb)
      common/ba/ y1,y2,y3,y4,y5,y6,y7

      real    ar(nxm,nzm),ai(nxm,nzm)
      common/rfft/ ar,ai

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer iv(nxi)
      real    s(nxq)
      integer i,i1,i2,id1,ie,ivi,ix,j,jj,jmm,jmmn,k,ks
      integer limit,lo2,lo4,mm,n2
      real    cosv,pint,sinv

      save

c
      if(iflg.eq.1) go to 100
      ie=0
      i=1
    6 ie=ie+1
      i=i+i
      if(i.lt.l) go to 6
      pint=8.*atan(1.)/float(l)
      limit=l/4+1
      do 7 i=1,limit
    7  s(i)=sin((i-1)*pint)
      iv(1)=0
      do 5 i=2,l
       ivi=0
       i1=i-1
       i2=i1/2
       do 8 j=1,ie
        ivi=ivi+ivi+i1-i2-i2
        i1=i2
        i2=i2/2
    8  continue
       iv(i)=ivi
    5 continue
      return
  100  lo2=l/2
       lo4=l/4
       n2=lo2
      do 1500 i=1,n2
       do 150 k=1,kl2
        y1(k)=ar(i,k)+ar(i+n2,k)
        y2(k)=ai(i,k)+ai(i+n2,k)
        y3(k)=ai(i,k)-ai(i+n2,k)
        y4(k)=ar(i,k)-ar(i+n2,k)
        ar(i,k)=y1(k)
        ai(i,k)=y2(k)
        ai(i+n2,k)=y3(k)
        ar(i+n2,k)=y4(k)
  150 continue
 1500 continue
       mm=l
       jj=1
  500  mm=mm/2
       n2=n2/2
       jj=jj+jj
       ks=-1
      do 3000 j=1,jj
       ks=ks+1
       ix=iv(ks+ks+1)
       jmm=mm*(j-1)
       jmmn=jmm+n2
       if(lo4-ix) 20,15,15
   15   sinv=s(ix+1)
        cosv=s(lo4-ix+1)
       go to 25
   20   sinv=s(lo2-ix+1)
        cosv=-s(ix-lo4+1)
   25  continue
       do 300 i=1,n2
        do 30 k=1,kl2
         y1(k)=ar(i+jmmn,k)*cosv-ai(i+jmmn,k)*sinv
         y2(k)=ar(i+jmmn,k)*sinv+ai(i+jmmn,k)*cosv
         y3(k)=ar(i+jmm,k)+y1(k)
         y4(k)=ai(i+jmm,k)+y2(k)
         y5(k)=ar(i+jmm,k)
         y6(k)=ai(i+jmm,k)
         ar(i+jmmn,k)=y5(k)-y1(k)
         ai(i+jmmn,k)=y6(k)-y2(k)
         ar(i+jmm,k)=y3(k)
   30    ai(i+jmm,k)=y4(k)
  300  continue
 3000 continue
      if(n2.ge.2) go to 500
      do 600 i=1,l
        id1=iv(i)+1
        if(i.ge.id1) go to 600
       do 60 k=1,kl2
         y1(k)=ar(id1,k)
         y2(k)=ai(id1,k)
         ar(id1,k)=ar(i,k)
         ai(id1,k)=ai(i,k)
         ar(i,k)=y1(k)
         ai(i,k)=y2(k)
   60  continue
  600 continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundy (x,x1)
      parameter (NX=514,NZ=43)
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    x(nx,nz),x1(nx,nz)
c     ******   periodic lateral boundary condition
      do 10 k=1,kmax
       x(1,k)=x(iles,k)
       x1(1,k)=x1(iles,k)
       x(imax,k)=x(2,k)
       x1(imax,k)=x1(2,k)
   10 continue
cc    ****   set boundary condition in z-direction
      do 20 i=1,imax
       x(i,1)=x(i,2)
       x1(i,1)=x1(i,2)
       x(i,kmax)=x(i,kles)
   20  x1(i,kmax)=x1(i,kles)
      return
      end
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
       pt(i)=pt(i)+afcp*y1(i)
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
         pt(i)=pt(i)-afcp*(psmlt(i)+pgmlt(i))
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
        pt(i)=pt(i)+afcp*y1(i)
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
            pt(i)=pt(i)+ascp*pint(i)
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
        pt(i)=pt(i)+avcp*cnd(i)+ascp*dep(i)
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
             PT(I)=PT(I)+AVCP*CND(I)+ASCP*DEP(I) ! cccshie by tao 5/3/01
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
         pt(i)=pt(i)+ascp*y1(i)
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
          pt(i)=pt(i)-avcp*ern(i)
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
          pt(i)=pt(i)-ascp*y1(i)
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
c
      subroutine sat_zero
c     ***************************************
      parameter (NX=514,NZ=43)
      parameter (nx16=16*nx)
      common/bsat/ xx1(nx,nz)
      common/bsat1/ xx2(nx,nz)
      common/badv1/ xx3(nx,nz)
      common/b1c/ qc(nx,nz)
      common/b1r/ qr(nx,nz)
      common/b1i/ qi(nx,nz)
      common/b1s/ qs(nx,nz)
      common/b1g/ qg(nx,nz)
      common/ba/ y1(nx16)
       CMIN=1.e-20
       do 100 k=1,nz
       do 100 i=1,nx
          xx1(i,k)=0.
          xx2(i,k)=0.
          xx3(i,k)=0.
        if (qc(i,k) .le. CMIN) qc(i,k)=0.0
        if (qr(i,k) .le. CMIN) qr(i,k)=0.0
        if (qi(i,k) .le. CMIN) qi(i,k)=0.0
        if (qs(i,k) .le. CMIN) qs(i,k)=0.0
        if (qg(i,k) .le. CMIN) qg(i,k)=0.0
  100  continue
       do 10 i=1,nx16
          y1(i)=0.
   10  continue
      return
      end
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine terv (irsg,rho,fv)
      implicit none
      integer nx,nz,nx10
      parameter (NX=514,NZ=43)
      parameter (nx10=10*nx)
      integer irsg
      real    fv(nz),rho(nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      real    qrn(nx,nz),qcg(nx,nz),qcs(nx,nz)
      common/b1r/ qrn
      common/b1g/ qcg
      common/b1s/ qcs
      real    ww1(nx,nz)
      COMMON/BTV/ WW1

      real    y1(nx),y2(nx),y3(nx),vr(nx),vs(nx),vg(nx),y4(nx10)
      common/ba/ y1,y2,y3,vr,vs,vg,y4

      real    ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq

      real    zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc
      common/rterv/ zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k
      real    vgcr,vscf

      save


c     ******************************************************************
      if(irsg.ne.0) go to 1
      do 10 k=2,kles
      do 10 i=2,iles
       ww1(i,k)=0.
         y1(i)=.5*rho(k)*(qrn(i,k)+qrn(i,k-1))
        if (y1(i) .gt. 1.e-15) then
           vs(i)=sqrt( y1(i) )
           vg(i)=sqrt( vs(i) )
          vr(i)=fv(k)*(vr0+vr1*vg(i)+vr2*vs(i)+vr3*vg(i)*vs(i))
         ww1(i,k)=max(vr(i), 0.E0)
        endif
   10 continue
      return
    1 if(irsg.ne.1) go to 2
      do 20 k=2,kles
       vscf=vsc*fv(k)
      do 20 i=2,iles
       ww1(i,k)=0.
         y1(i)=.5*rho(k)*(qcs(i,k)+qcs(i,k-1))
        if (y1(i) .gt. 1.e-15) then
         vs(i)=vscf*y1(i)**bsq
         ww1(i,k)=max(vs(i), 0.E0)
        endif
   20 continue
      return
    2 do 30 k=2,kles
       vgcr=vgc*fv(k)
      do 30 i=2,iles
       ww1(i,k)=0.
         y1(i)=.5*rho(k)*(qcg(i,k)+qcg(i,k-1))
        if (y1(i) .gt. 1.e-15) then
          vg(i)=vgcr*y1(i)**bgq
         ww1(i,k)=max(vg(i), 0.E0)
        endif
   30 continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SATICEL (n,id,fv)
C     (LIN)  COMPUTE ICE PHASE MICROPHYSICS AND SATURATION PROCESSES
      PARAMETER (NX=514,NZ=43,NT=2880,ITT=244)
      PARAMETER (NT5=5*NT,NB=NX*NZ-29*NX,NB1=NX*NZ-13*NX)

      integer n,id
      real    fv(nz)

      common/timestat/ ndt_stat,itime_ave,mmave
      integer itoga,ISFC,ICE,ICE2
      common/itoga/ itoga,ISFC,ICE,ICE2
      common/iceopt/ ice913,ilif
      common/ilin/ lin
      common/iice/ new_ice_sat
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar

      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bsize/ it1,it2,kt1,kt2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141

      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC
      COMMON/BSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,BND1,RN2,BND2,
     1   RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,RN101,RN10A,RN11,RN11A,
     2   RN12,RN12A(31),RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN17,
     3   RN17A,RN17B,RN17C,RN18,RN18A,RN19,RN19A,RN19B,RN20,RN20A,RN20B,
     4   BND3,RN21,RN22,RN23,RN23A,RN23B,RN25,RN25A(31),RN30A,RN30B,
     5   RN30C,RN31,BETA,RN32

      COMMON/B1T/ DPT(NX,NZ)
      COMMON/B1Q/ DQV(NX,NZ)
      COMMON/B1C/ QCL(NX,NZ)
      COMMON/B1R/ QRN(NX,NZ)
      COMMON/B2T/ DPT1(NX,NZ)
      COMMON/B2Q/ DQV1(NX,NZ)
      COMMON/B2C/ QCL1(NX,NZ)
      COMMON/B2R/ QRN1(NX,NZ)
      COMMON/B1I/ QCI(NX,NZ)
      COMMON/B1S/ QCS(NX,NZ)
      COMMON/B1G/ QCG(NX,NZ)
      COMMON/B2I/ QCI1(NX,NZ)
      COMMON/B2S/ QCS1(NX,NZ)
      COMMON/B2G/ QCG1(NX,NZ)
      COMMON/B2W/ WW1(NX,NZ)
      COMMON/SLWAVE/ RSW(NX,NZ),RLW(NX,NZ)

      COMMON/BSAT/ WGACR(NX),DEP(NX),RGMP(NX),
     1   DD(NX),DD1(NX),QVS(NX),DM(NX),RSUB1(NX),COL(NX),CND(NX),RQ(NX),
     2   ERN(NX),SCV(NX),TCA(NX),DWV(NX),ZR(NX),VR(NX),ZS(NX),VS(NX),
     3   ZG(NX),VG(NX),EGS(NX),ESI(NX),QSI(NX),SSI(NX),QSW(NX),SSW(NX),
     4   PIHOM(NX),PIDW(NX),DDC(NB)

      COMMON/BSAT1/ QSACW(NX),PRACI(NX),PIACR(NX),PRAUT(NX),
     1   PRACW(NX),PSFW(NX),PSFI(NX),DGACS(NX),DGACW(NX),DGACI(NX),
     2   DGACR(NX),PGACS(NX),WGACS(NX),QGACW(NX),WGACI(NX),QGACR(NX),
     3   PGWET(NX),PGAUT(NX),PRACS(NX),PSACR(NX),QSACR(NX),PGFR(NX),
     4   PSMLT(NX),PGMLT(NX),PSDEP(NX),PSSUB(NX),PGSUB(NX),PINT(NX),
     5   PIDEP(NX),DDA(NB)


      COMMON/BADV/ PT(NX),QV(NX),QC(NX),QR(NX),QI(NX),QS(NX),QG(NX),
     1   XX0(NX),XX00(NX),PIMLT(NX),PSAUT(NX),PSACI(NX),PSACW(NX),
     2   DDB(NB1)

      real    y1(nx),y2(nx),y3(nx),y4(nx),y5(nx),tair(nx),tairc(nx),
     $   pr(nx),ps(nx),pg(nx),prn(nx),psn(nx),dlt1(nx),dlt2(nx),
     $   dlt3(nx),rtair(nx)
      common/ba/ y1,y2,y3,y4,y5,tair,tairc,pr,ps,pg,prn,psn,dlt1,dlt2,
     $           dlt3,rtair
c
c
      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     $   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    srro(nz),qrro(nz),sqc(nz),sqr(nz),sqi(nz),sqs(nz),
     $   sqg(nz),stqc(nz),stqr(nz),stqi(nz),stqs(nz),stqg(nz),ttt(nt5)
      common/brh1/ srro,qrro,sqc,sqr,sqi,sqs,
     $             sqg,stqc,stqr,stqi,stqs,stqg,ttt

      COMMON/B7/ SQAAQ(NZ),SQTT(NZ),SQVV(NZ),SQAAT(NZ),SQAK(NZ)
      real    y0(nx),ths(nx,itt),qsTT(nx,itt),tsTT(nx,itt),pss(nx,itt)
      common/bls/ y0,ths,qsTT,tsTT,pss

      integer it(nx),ics(nx,4),ibz(nx,4),iv(nt)
      common/bch/ it,iv,ics,ibz
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON/BSTS/ THOM(NZ,4,7),TDW(NZ,4,7),TMLT(NZ,4,7),
     1 SAUT(NZ,4,7),SACI(NZ,4,7),SACW(NZ,4,7),RACI(NZ,4,7),
     2 TACR(NZ,4,7),RAUT(NZ,4,7),RACW(NZ,4,7),SFW(NZ,4,7),
     3 SFI(NZ,4,7),GACS(NZ,4,7),GACW(NZ,4,7),
     4 GACI(NZ,4,7),GACR(NZ,4,7),GWET(NZ,4,7),GAUT(NZ,4,7),
     5 RACS(NZ,4,7),SACR(NZ,4,7),GFR(NZ,4,7),SMLT(NZ,4,7),
     6 GMLT(NZ,4,7),SDEP(NZ,4,7),SSUB(NZ,4,7),GSUB(NZ,4,7),
     7 PERN(NZ,4,7),D3RI(NZ,4,7),D3IR(NZ,4,7),D2SR(NZ,4,7),
     8 D2RS(NZ,4,7),GDRY(NZ,4,7),COC(NZ,4,7),COE(NZ,4,7),
     9 SMF0(NZ,4,7),QC0(NZ,4,7),QR0(NZ,4,7),QI0(NZ,4,7),
     1 QS0(NZ,4,7),QG0(NZ,4,7),SQC0(NZ,4,7),SQR0(NZ,4,7),
     2 SQI0(NZ,4,7),SQS0(NZ,4,7),SQG0(NZ,4,7),ERNS(NZ,4,7),
     3 WGRS(NZ,4,7),QSWS(NZ,4,7),TBT(NZ,4),QBT(NZ,4)


      COMMON/BSTS4/ SRSW(NZ,4,7),SRLW(NZ,4,7)

      COMMON/BSTSLIN/ qv0l(nz,4,7),tt0l(nz,4,7),sgpt(nz,4,7),
     1  tsqq(nz,4,7),stt0l(nz,4,7),sqv0l(nz,4,7),ssgpt(nz,4,7),
     2  tsqq1(nz,4,7),q1tl(nz,4,7),sqqdtl(nz,4,7),sqdtl(nz,4,7),
     3  stdtl(nz,4,7)

      real   S9(16,NZ),S10(16,NZ),S11(16,NZ),S12(16,NZ),S13(16,NZ),
     $       S14(16,NZ),S15(16,NZ),S16(16,NZ),S17(16,NZ),S18(16,NZ),
     $       S19(16,NZ),S20(16,NZ),S21(16,NZ),SCNT(16,NZ),SN9(5,NZ),
     $       SN10(5,NZ),SN11(5,NZ),SN12(5,NZ),SN13(5,NZ),SN14(5,NZ),
     $       SN15(5,NZ),SN16(5,NZ),SN17(5,NZ),SN18(5,NZ),SN19(5,NZ),
     $       SN20(5,NZ),SN21(5,NZ),SNCNT(5,NZ),SCU1(NZ),SED1(NZ)
      COMMON/BCS/ S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,
     $            S19,S20,S21,SCNT,SN9,SN10,SN11,SN12,SN13,SN14,
     $            SN15,SN16,SN17,SN18,SN19,SN20,SN21,SNCNT,SCU1,SED1

      real    SS9(16,NZ),SS10(16,NZ),SS11(16,NZ),SS12(16,NZ),
     $ SS13(16,NZ),SS14(16,NZ),SS15(16,NZ),SS16(16,NZ),SS17(16,NZ),
     $ SS18(16,NZ),SS19(16,NZ),SS20(16,NZ),SS21(16,NZ),SSCNT(16,NZ),
     $ SSN9(5,NZ),SSN10(5,NZ),SSN11(5,NZ),SSN12(5,NZ),SSN13(5,NZ),
     $ SSN14(5,NZ),SSN15(5,NZ),SSN16(5,NZ),SSN17(5,NZ),SSN18(5,NZ),
     $ SSN19(5,NZ),SSN20(5,NZ),SSN21(5,NZ),SSNCNT(5,NZ)
      COMMON/BCSS/ SS9,SS10,SS11,SS12,SS13,SS14,SS15,SS16,SS17,
     $         SS18,SS19,SS20,SS21,SSCNT,SSN9,SSN10,SSN11,SSN12,SSN13,
     $         SSN14,SSN15,SSN16,SSN17,SSN18,SSN19,SSN20,SSN21,SSNCNT
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    rby(7)
      common/bch1/ rby

      DIMENSION AA1(31),AA2(31),dbz(nx)
      DATA AA1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,
     1   .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,
     2   .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,
     3   .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,
     4   .6513E-6,.5956E-6,.5333E-6,.4834E-6/
      DATA AA2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
      save
C     ******   THREE CLASSES OF ICE-PHASE   (LIN ET AL, 1983)  *********
       RT0=1./(T0-T00)
       D22T=D2T
       IF(IJKADV .EQ. 0) THEN
         D2T=D2T
       ELSE
         D2T=DT
       ENDIF
c
        UCOR=3071.29/TNW**.75                                            
        UCOS=687.97*ROQS**.25/TNS**.75
        UCOG=687.97*ROQG**.25/TNG**.75
        UWET=4.464**.95
       cmin=1.e-20
       cmin1=1.e-12
c
       A0=ndt_stat*.5*RIL2
      DO 10 I=1,imax
   10  IT(I)=1
       FT=DT/D2T
       RFT=RIL2*FT
        BW3=BW+3.
        BS3=BS+3.
        BG3=BG+3.
        BWH5=2.5+BWH
        BSH5=2.5+BSH
        BGH5=2.5+BGH
        BW6=BW+6.
        BETAH=.5*BETA
       R10T=RN10*D2T
       R11AT=RN11A*D2T
       R19BT=RN19B*D2T
       R20T=-RN20*D2T
       R23T=-RN23*D2T
       R25A=RN25
CC    ******************************************************************
      DO 1000 K=2,KLES
c        if (ijkadv .eq. 1) then
c          tb0=ta(k)
c          qb0=qa(k)
c        else
          tb0=ta1(k)
          qb0=qa1(k)
c        endif
c       P00=P0(K)
       RP0=3.799052E3/P0(K)
       PI0=PI(K)
       PIR=1./(PI(K))
       PR0=1./P0(K)
       R00=RHO(K)
       RR0=RRHO(K)
       RRS=SRRO(K)
       RRQ=QRRO(K)
       FV0=FV(K)
       FVS=SQRT(FV(K))
         ZRR=1.E5*ZRC*RRQ
         ZSR=1.E5*ZSC*RRQ
         ZGR=1.E5*ZGC*RRQ
         CP409=C409*PI0
         CV409=C409*AVC
         CP580=C580*PI0
         CS580=C580*ASC
         ALVR=R00*ALV
         AFCP=AFC*PIR
         AVCP=AVC*PIR
         ASCP=ASC*PIR
         VRCF=VRC*FV0
         VSCF=VSC*FV0

         VGCR=VGC*RRS
         VGCF=VGC*FV0

         DWVP=C879*PR0
         R3F=RN3*FV0
         R4F=RN4*FV0
         R5F=RN5*FV0
         R6F=RN6*FV0
         R7R=RN7*RR0
         R8R=RN8*RR0
         R9R=RN9*RR0
         R101F=RN101*FVS
         R10AR=RN10A*R00
         R11RT=RN11*RR0*D2T
         R12R=RN12*R00
         R14R=RN14*RRS
         R15R=RN15*RRS
         R15AR=RN15A*RRS
         R16R=RN16*RR0
         R17R=RN17*RR0
         R17AQ=RN17A*RRQ
         R18R=RN18*RR0
         R19RT=RN19*RR0*D2T
         R19AQ=RN19A*RRQ
         R20BQ=RN20B*RRQ
         R22F=RN22*FV0
         R23AF=RN23A*FVS
         R23BR=RN23B*R00
         R25RT=RN25*RR0*D2T
         R31R=RN31*RR0
         R32RT=RN32*D2T*RRS
       DO 100 I=2,ILES
        PT(I)=DPT(I,K)
        QV(I)=DQV(I,K)
        QC(I)=QCL(I,K)
        QR(I)=QRN(I,K)
        QI(I)=QCI(I,K)
        QS(I)=QCS(I,K)
        QG(I)=QCG(I,K)
c        IF (QV(I)+QB0 .LE. 0.) QV(I)=-QB0
        IF (QC(I) .LE. cmin) QC(I)=0.
        IF (QR(I) .LE. cmin) QR(I)=0.
        IF (QI(I) .LE. cmin) QI(I)=0.
        IF (QS(I) .LE. cmin) QS(I)=0.
        IF (QG(I) .LE. cmin) QG(I)=0.
        TAIR(I)=(PT(I)+TB0)*PI0
        TAIRC(I)=TAIR(I)-T0
C     ***   COMPUTE ZR,ZS,ZG,VR,VS,VG      *****************************
          ZR(I)=ZRR
          ZS(I)=ZSR
          ZG(I)=ZGR
          VR(I)=0.0
          VS(I)=0.0
          VG(I)=0.0
           IF (QR(I) .GT. cmin) THEN
             DD(I)=R00*QR(I)
              Y1(I)=DD(I)**.25
            ZR(I)=ZRC/Y1(I)
            VR(I)=min(max(VRCF*DD(I)**BWQ, 0.0), 1500.)
           ENDIF
           IF (QS(I) .GT. cmin) THEN
             DD(I)=R00*QS(I)
              Y1(I)=DD(I)**.25
            ZS(I)=ZSC/Y1(I)
            VS(I)=min(max(VSCF*DD(I)**BSQ, 0.0), 400.)
           ENDIF
           IF (QG(I) .GT. cmin) THEN
             DD(I)=R00*QG(I)
              Y1(I)=DD(I)**.25
            ZG(I)=ZGC/Y1(I)
            IF (LIN .EQ. 1) THEN
             vg(i)=min(max(vgcr*dd(i)**bgq, 0.), 1500.)
            ELSE
             vg(i)=min(max(vgcf*dd(i)**bgq, 0.), 800.)
            ENDIF
           ENDIF
         IF (QR(I) .LE. cmin1) VR(I)=0.0
         IF (QS(I) .LE. cmin1) VS(I)=0.0
         IF (QG(I) .LE. cmin1) VG(I)=0.0
C     ******************************************************************
C     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
C     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
C     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
C     ***   Y2 : KINETIC VISCOSITY (V)
          Y1(I)=C149*TAIR(I)**1.5/(TAIR(I)+120.)
        DWV(I)=DWVP*TAIR(I)**1.81
        TCA(I)=C141*Y1(I)
        SCV(I)=1./((RR0*Y1(I))**.1666667*DWV(I)**.3333333)
  100  CONTINUE
C*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
C*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
C*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
C*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
C*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
        DO 125 I=2,ILES
         PSAUT(I)=0.0
         PSACI(I)=0.0
         PRACI(I)=0.0
         PIACR(I)=0.0
         PSACW(I)=0.0
         QSACW(I)=0.0
          DD(I)=1./ZS(I)**BS3
         IF (TAIR(I).LT.T0) THEN
           ESI(I)=EXP(.025*TAIRC(I))
          PSAUT(I)=MAX(RN1*ESI(I)*(QI(I)-BND1) ,0.0)
          PSACI(I)=R3F*ESI(I)*QI(I)*DD(I)
          PSACW(I)=R4F*QC(I)*DD(I)
          PRACI(I)=R5F*QI(I)/ZR(I)**BW3
          PIACR(I)=R6F*QI(I)/ZR(I)**BW6
         ELSE
          QSACW(I)=R4F*QC(I)*DD(I)
         ENDIF
C* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
C* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
CM        PRAUT(I)=MAX((RD1*(QC(I)-BOUND)), 0.0)
CM        PRACW(I)=0.0
CM         IF (QR(I).GT.1.E-40) THEN
CM          PRACW(I)=RD2*QC(I)*QR(I)**.875
CM         ENDIF
          PRACW(I)=R22F*QC(I)/ZR(I)**BW3
          PRAUT(I)=0.0
            Y1(I)=QC(I)-BND3
           IF (Y1(I).GT.0.0) THEN
            PRAUT(I)=R00*Y1(I)*Y1(I)/(1.2E-4+RN21/Y1(I))
           ENDIF
C* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
C* 13 * PSFI : BERGERON PROCESSES FOR QS                          **13**
          PSFW(I)=0.0
          PSFI(I)=0.0
           IF(TAIR(I).LT.T0) THEN
             Y1(I)=MAX( MIN(TAIRC(I), -1.), -31.)
             IT(I)=INT(ABS(Y1(I)))
             Y1(I)=RN12A(IT(I))
             Y2(I)=RN12B(IT(I))
            PSFW(I)=MAX(D2T*Y1(I)*(Y2(I)+R12R*QC(I))*QI(I),0.0)
            PSFI(I)=RN13(IT(I))*QI(I)
           ENDIF
CTTT***** QG=QG+MIN(PGDRY,PGWET)
C*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
C* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
C* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
           EE1=1.
           EE2=0.09
          EGS(I)=EE1*EXP(EE2*TAIRC(I))
           IF (TAIR(I).GE.T0) EGS(I)=1.0
          Y1(I)=ABS(VG(I)-VS(I))
          Y2(I)=ZS(I)*ZG(I)
          Y3(I)=5./Y2(I)
          Y4(I)=.08*Y3(I)*Y3(I)
          Y5(I)=.05*Y3(I)*Y4(I)
          DD(I)=Y1(I)*(Y3(I)/ZS(I)**5+Y4(I)/ZS(I)**3+Y5(I)/ZS(I))
         PGACS(I)=R9R*EGS(I)*DD(I)
         DGACS(I)=PGACS(I)
         WGACS(I)=R9R*DD(I)
          Y1(I)=1./ZG(I)**BG3
         DGACW(I)=MAX(R14R*QC(I)*Y1(I), 0.0)
         QGACW(I)=DGACW(I)
          Y1(I)=ABS(VG(I)-VR(I))
          Y2(I)=ZR(I)*ZG(I)
          Y3(I)=5./Y2(I)
          Y4(I)=.08*Y3(I)*Y3(I)
          Y5(I)=.05*Y3(I)*Y4(I)
           DD(I)=R16R*Y1(I)*(Y3(I)/ZR(I)**5+Y4(I)/ZR(I)**3+Y5(I)/ZR(I))
         DGACR(I)=MAX(DD(I), 0.0)
         QGACR(I)=DGACR(I)
         IF (TAIR(I).GE.T0) THEN
          DGACS(I)=0.0
          WGACS(I)=0.0
          DGACW(I)=0.0
          DGACR(I)=0.0
         ELSE
          PGACS(I)=0.0
          QGACW(I)=0.0
          QGACR(I)=0.0
         ENDIF
C*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
C* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
C* 17 * PGWET : WET GROWTH OF QG                                  **17**
          DGACI(I)=0.0
          WGACI(I)=0.0
          PGWET(I)=0.0
          IF (TAIR(I).LT.T0 .AND. TAIR(IJ) .GE. 253.16) THEN
             Y1(I)=QI(I)/ZG(I)**BG3
            DGACI(I)=R15R*Y1(I)
            WGACI(I)=R15AR*Y1(I)
             Y1(I)=1./(ALF+RN17C*TAIRC(I))
             Y3(I)=.78/ZG(I)**2+R17AQ*SCV(I)/ZG(I)**BGH5
             Y4(I)=ALVR*DWV(I)*(RP0-(QV(I)+QB0))-TCA(I)*TAIRC(I)
            DD(I)=Y1(I)*(R17R*Y4(I)*Y3(I)
     1                      +(WGACI(I)+WGACS(I))*(ALF+RN17B*TAIRC(I)))
           PGWET(I)=MAX(DD(I), 0.0)
          ENDIF
  125   CONTINUE
C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
C********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
        DO 150 I=2,ILES

           Y1(I)=QC(I)/D2T
          PSACW(I)=MIN(Y1(I), PSACW(I))
          PRAUT(I)=MIN(Y1(I), PRAUT(I))
          PRACW(I)=MIN(Y1(I), PRACW(I))
          PSFW(I)= MIN(Y1(I), PSFW(I))
          DGACW(I)=MIN(Y1(I), DGACW(I))
          QSACW(I)=MIN(Y1(I), QSACW(I))
          QGACW(I)=MIN(Y1(I), QGACW(I))

        Y1(I)=(PSACW(I)+PRAUT(I)+PRACW(I)+PSFW(I)+DGACW(I)+QSACW(I)
     1         +QGACW(I))*D2T
       QC(I)=QC(I)-Y1(I)
        if (qc(i).lt.0.0) then
          a1=1.
           if(y1(i) .ne. 0.) a1=qc(i)/y1(i)+1.
         psacw(i)=psacw(i)*a1
         praut(i)=praut(i)*a1
         pracw(i)=pracw(i)*a1
         psfw(i)=psfw(i)*a1
         dgacw(i)=dgacw(i)*a1
         qsacw(i)=qsacw(i)*a1
         qgacw(i)=qgacw(i)*a1
         qc(i)=0.0
        endif



            Y1(I)=QI(I)/D2T
           PSAUT(I)=MIN(Y1(I), PSAUT(I))
           PSACI(I)=MIN(Y1(I), PSACI(I))
           PRACI(I)=MIN(Y1(I), PRACI(I))
           PSFI(I)= MIN(Y1(I), PSFI(I))
           DGACI(I)=MIN(Y1(I), DGACI(I))
           WGACI(I)=MIN(Y1(I), WGACI(I))

        Y1(I)=(PSAUT(I)+PSACI(I)+PRACI(I)+PSFI(I)+DGACI(I)+WGACI(I))*D2T
        QI(I)=QI(I)-Y1(I)

        if (qi(i).lt.0.0) then
          a2=1.
           if(y1(i) .ne. 0.) a2=qi(i)/y1(i)+1.
         psaut(i)=psaut(i)*a2
         psaci(i)=psaci(i)*a2
         praci(i)=praci(i)*a2
         psfi(i)=psfi(i)*a2
         dgaci(i)=dgaci(i)*a2
         wgaci(i)=wgaci(i)*a2
         qi(i)=0.0
        endif
c
C******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
         WGACR(I)=PGWET(I)-DGACW(I)-WGACI(I)-WGACS(I)
          Y2(I)=DGACW(I)+DGACI(I)+DGACR(I)+DGACS(I)
         IF (PGWET(I).GE.Y2(I)) THEN
          WGACR(I)=0.0
          WGACI(I)=0.0
          WGACS(I)=0.0
         ELSE
          DGACR(I)=0.0
          DGACI(I)=0.0
          DGACS(I)=0.0
         ENDIF
c******************************************************
c
         DLT3(I)=0.0
         DLT2(I)=0.0
          IF (TAIR(I).LT.T0) THEN
           IF (QR(I).LT.1.E-4) THEN
             DLT3(I)=1.0
             DLT2(I)=1.0
           ENDIF
           IF (QS(I).GE.1.E-4) THEN
             DLT2(I)=0.0
           ENDIF
          ENDIF
         PR(I)=(QSACW(I)+PRAUT(I)+PRACW(I)+QGACW(I))*D2T
         PS(I)=(PSAUT(I)+PSACI(I)+PSACW(I)+PSFW(I)+PSFI(I)
     1          +DLT3(I)*PRACI(I))*D2T
         PG(I)=((1.-DLT3(I))*PRACI(I)+DGACI(I)+WGACI(I)+DGACW(I))*D2T
  150   CONTINUE
C*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
C*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
        DO 175 I=2,ILES
          Y1(I)=ABS(VR(I)-VS(I))
          Y2(I)=ZR(I)*ZS(I)
          Y3(I)=5./Y2(I)
          Y4(I)=.08*Y3(I)*Y3(I)
          Y5(I)=.05*Y3(I)*Y4(I)
         PRACS(I)=R7R*Y1(I)*(Y3(I)/ZS(I)**5+Y4(I)/ZS(I)**3+Y5(I)/ZS(I))
         PSACR(I)=R8R*Y1(I)*(Y3(I)/ZR(I)**5+Y4(I)/ZR(I)**3+Y5(I)/ZR(I))
         QSACR(I)=PSACR(I)
         IF (TAIR(I).GE.T0) THEN
          PRACS(I)=0.0
          PSACR(I)=0.0
         ELSE
          QSACR(I)=0.0
         ENDIF
C*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
C* 18 * PGFR : FREEZING OF QR TO QG                               **18**
         PGAUT(I)=0.0
         PGFR(I)=0.0
          IF (TAIR(I) .LT. T0) THEN
            Y1(I)=EXP(.09*TAIRC(I))
            IF (LIN .EQ. 1) PGAUT(IJ)=MAX(RN2*Y1(IJ)*(QS(IJ)-BND2), 0.0)
            Y2(I)=EXP(RN18A*(T0-TAIR(I)))
           PGFR(I)=MAX(R18R*(Y2(I)-1.)/ZR(I)**7., 0.0)
          ENDIF
  175   CONTINUE
C********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
C********   HANDLING THE NEGATIVE SNOW (QS)          *******************
        DO 200 I=2,ILES

          Y1(I)=QR(I)/D2T
         PIACR(I)=MIN(Y1(I), PIACR(I))
         DGACR(I)=MIN(Y1(I), DGACR(I))
         WGACR(I)=MIN(Y1(I), WGACR(I))
         PSACR(I)=MIN(Y1(I), PSACR(I))
         PGFR(I)= MIN(Y1(I), PGFR(I))

          Y1(I)=(PIACR(I)+DGACR(I)+WGACR(I)+PSACR(I)+PGFR(I))*D2T
         QR(I)=QR(I)+PR(I)-Y1(I)
         if (qr(i).lt.0.0) then
           a1=1.
            if(y1(i) .ne. 0.) a1=qr(i)/y1(i)+1.
          PIACR(I)=PIACR(I)*A1
          DGACR(I)=DGACR(I)*A1
          WGACR(I)=WGACR(I)*A1
          PGFR(I)=PGFR(I)*A1
          PSACR(I)=PSACR(I)*A1
          QR(I)=0.0
         ENDIF
          PRN(I)=D2T*((1.-DLT3(I))*PIACR(I)+DGACR(I)+WGACR(I)
     1                +(1.-DLT2(I))*PSACR(I)+PGFR(I))
          PS(I)=PS(I)+D2T*(DLT3(I)*PIACR(I)+DLT2(I)*PSACR(I))
           PRACS(I)=(1.-DLT2(I))*PRACS(I)

           Y1(I)=QS(I)/D2T
          PGACS(I)=MIN(Y1(I), PGACS(I))
          DGACS(I)=MIN(Y1(I), DGACS(I))
          WGACS(I)=MIN(Y1(I), WGACS(I))
          PGAUT(I)=MIN(Y1(I), PGAUT(I))
          PRACS(I)=MIN(Y1(I), PRACS(I))

          PSN(I)=D2T*(PGACS(I)+DGACS(I)+WGACS(I)+PGAUT(I)+PRACS(I))
         QS(I)=QS(I)+PS(I)-PSN(I)
          if (qs(i).lt.0.) then
            A2=1.
             IF (PSN(I) .NE. 0.0) a2=qs(i)/psn(i)+1.
           pgacs(i)=pgacs(i)*a2
           dgacs(i)=dgacs(i)*a2
           wgacs(i)=wgacs(i)*a2
           pgaut(i)=pgaut(i)*a2
           pracs(i)=pracs(i)*a2
c           psn(i)=psn(i)*a2
           qs(i)=0.0
          endif
c
          PSN(I)=D2T*(PGACS(I)+DGACS(I)+WGACS(I)+PGAUT(I)+PRACS(I))
           Y2(I)=D2T*(PSACW(I)+PSFW(I)+DGACW(I)+PIACR(I)+DGACR(I)
     1                +WGACR(I)+PSACR(I)+PGFR(I))
         PT(I)=PT(I)+AFCP*Y2(I)
         QG(I)=QG(I)+PG(I)+PRN(I)+PSN(I)
  200   CONTINUE
C* 11 * PSMLT : MELTING OF QS                                     **11**
C* 19 * PGMLT : MELTING OF QG TO QR                               **19**
        DO 225 I=2,ILES
        PSMLT(I)=0.0
        PGMLT(I)=0.0
           TAIR(I)=(PT(I)+TB0)*PI0
        IF (TAIR(I).GE.T0) THEN
           TAIRC(I)=TAIR(I)-T0
           Y1(I)=TCA(I)*TAIRC(I)-ALVR*DWV(I)*(RP0-(QV(I)+QB0))
           Y2(I)=.78/ZS(I)**2+R101F*SCV(I)/ZS(I)**BSH5
           DD(I)=R11RT*Y1(I)*Y2(I)+R11AT*TAIRC(I)*(QSACW(I)+QSACR(I))
         PSMLT(I)=MAX(0.0, MIN(DD(I), QS(I)))
           Y3(I)=.78/ZG(I)**2+R19AQ*SCV(I)/ZG(I)**BGH5
           DD1(I)=R19RT*Y1(I)*Y3(I)+R19BT*TAIRC(I)*(QGACW(I)+QGACR(I))
         PGMLT(I)=MAX(0.0, MIN(DD1(I), QG(I)))
C         IF(ICS5(I,3,1).EQ.1 .OR. ICS5(I,4,1).EQ.1) THEN
C         IF(ICS5(I,2,1).EQ.1) THEN
C           PSMLT(I)=0.0
C           PGMLT(I)=0.0
C         ENDIF
         PT(I)=PT(I)-AFCP*(PSMLT(I)+PGMLT(I))
         QR(I)=QR(I)+PSMLT(I)+PGMLT(I)
         QS(I)=QS(I)-PSMLT(I)
         QG(I)=QG(I)-PGMLT(I)
        ENDIF

C* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
C* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0 < T <= T00)     **25**
C* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
        IF (QC(I).LE.cmin) QC(I)=0.0
        IF (QI(I).LE.cmin) QI(I)=0.0
          TAIR(I)=(PT(I)+TB0)*PI0
         PIHOM(I)=CVMGP(0.,QC(I),TAIR(I)-T00)
         PIMLT(I)=CVMGP(QI(I),0.,TAIR(I)-T0)
         PIDW(I)=0.0
         IF (TAIR(I).LT.T0 .AND. TAIR(I).GT.T00) THEN
          TAIRC(I)=TAIR(I)-T0
            Y1(I)=MAX( MIN(TAIRC(I), -1.), -31.)
            IT(I)=INT(ABS(Y1(I)))
CLIN        Y2(I)=RN25A(IT(I))
CLIN        Y3(I)=EXP(.5*ABS(TAIRC(I)))
CLIN       PIDW(I)=MIN(R25RT*Y2(I)*Y3(I), QC(I))
C
           Y2(I)=AA1(IT(I))
           Y3(I)=AA2(IT(I))
           Y4(I)=EXP(.5*ABS(TAIRC(I)))
           DD(I)=(R00*QI(I)/(R25A*Y4(I)))**Y3(I)
          PIDW(I)=MIN(R25RT*Y2(I)*Y4(I)*DD(I), QC(I))
         ENDIF
          Y1(I)=PIHOM(I)-PIMLT(I)+PIDW(I)
        PT(I)=PT(I)+AFCP*Y1(I)
        QC(I)=QC(I)-Y1(I)
        QI(I)=QI(I)+Y1(I)
C* 31 * PINT  : INITIATION OF QI                                  **31**
C* 32 * PIDEP : DEPOSITION OF QI                                  **32**
        PINT(I)=0.0
         TAIR(I)=(PT(I)+TB0)*PI0
         IF (TAIR(I) .LT. T0) THEN
            if (qi(i) .le. cmin) qi(i)=0.
           TAIRC(I)=TAIR(I)-T0
           DD(I)=R31R*EXP(BETA*TAIRC(I))
            RTAIR(I)=1./(TAIR(I)-C76)
             Y2(I)=EXP(C218-C580*RTAIR(I))
            QSI(I)=RP0*Y2(I)
            ESI(I)=C610*Y2(I)
            SSI(I)=(QV(I)+QB0)/QSI(I)-1.
           DM(I)=MAX( (QV(I)+QB0-QSI(I)), 0.)
            RSUB1(I)=CS580*QSI(I)*RTAIR(I)*RTAIR(I)
          PINT(I)=max(MIN(DD(I),DM(I)), 0.)
             Y1(I)=1./TAIR(I)
             Y2(I)=EXP(BETAH*TAIRC(I))
             Y3(I)=SQRT(QI(I))
           DD(I)=Y1(I)*(RN30A*Y1(I)-RN30B)+RN30C*TAIR(I)/ESI(I)
          PIDEP(I)=MAX(R32RT*SSI(I)*Y2(I)*Y3(I)/DD(I), 0.0)
          PINT(I)=PINT(I)+PIDEP(I)
           DEP(I)=DM(I)/(1.+RSUB1(I))
           PINT(I)=MIN(PINT(I), DEP(I))
cc            if (PINT(i) .le. cmin) PINT(i)=0.
          PT(I)=PT(I)+ASCP*PINT(I)
          QV(I)=QV(I)-PINT(I)
          QI(I)=QI(I)+PINT(I)
         ENDIF
  225   CONTINUE
C*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
        DO 250 I=2,ILES
         TAIR(I)=(PT(I)+TB0)*PI0
        CND(I)=RT0*(TAIR(I)-T00)
        DEP(I)=RT0*(T0-TAIR(I))
          Y1(I)=1./(TAIR(I)-C358)
          Y2(I)=1./(TAIR(I)-C76)
         QSW(I)=RP0*EXP(C172-C409*Y1(I))
         QSI(I)=RP0*EXP(C218-C580*Y2(I))
          DD(I)=CP409*Y1(I)*Y1(I)
          DD1(I)=CP580*Y2(I)*Y2(I)
         IF (QC(I).LE.cmin) QC(I)=cmin
         IF (QI(I).LE.cmin) QI(I)=cmin
         IF (TAIR(I).GE.T0) THEN
          DEP(I)=0.0
          CND(I)=1.0
          QI(I)=0.0
         ENDIF
         IF (TAIR(I).LE.T00) THEN
          CND(I)=0.0
          DEP(I)=1.0
          QC(I)=0.0
         ENDIF
          Y5(I)=AVCP*CND(I)+ASCP*DEP(I)
           Y1(I)=QC(I)*QSW(I)/(QC(I)+QI(I))
           Y2(I)=QI(I)*QSI(I)/(QC(I)+QI(I))
          Y4(I)=DD(I)*Y1(I)+DD1(I)*Y2(I)
         QVS(I)=Y1(I)+Y2(I)
         DM(I)=QV(I)+QB0-QVS(I)
         RSUB1(I)=DM(I)/(1.+Y4(I)*Y5(I))
        CND(I)=CND(I)*RSUB1(I)
        DEP(I)=DEP(I)*RSUB1(I)
         IF (QC(I).LE.cmin) QC(I)=0.
         IF (QI(I).LE.cmin) QI(I)=0.
CC    ******   CONDENSATION OR EVAPORATION OF QC  ******
         CND(I)=MAX(-QC(I), CND(I))
CC    ******   DEPOSITION OR SUBLIMATION OF QI    ******
         DEP(I)=MAX(-QI(I), DEP(I))
        PT(I)=PT(I)+AVCP*CND(I)+ASCP*DEP(I)
        QV(I)=QV(I)-CND(I)-DEP(I)
        QC(I)=QC(I)+CND(I)
        QI(I)=QI(I)+DEP(I)
  250   CONTINUE
C* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
C* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
        DO 275 I=2,ILES
         PSDEP(I)=0.0
         PSSUB(I)=0.0
         PGSUB(I)=0.0
          TAIR(I)=(PT(I)+TB0)*PI0
          IF(TAIR(I).LT.T0) THEN
             IF(QS(I).LT.cmin) QS(I)=0.0
             IF(QG(I).LT.cmin) QG(I)=0.0
             RTAIR(I)=1./(TAIR(I)-C76)
            QSI(I)=RP0*EXP(C218-C580*RTAIR(I))
            SSI(I)=(QV(I)+QB0)/QSI(I)-1.
              Y1(I)=R10AR/(TCA(I)*TAIR(I)**2)+1./(DWV(I)*QSI(I))
              Y2(I)=.78/ZS(I)**2+R101F*SCV(I)/ZS(I)**BSH5
           PSDEP(I)=R10T*SSI(I)*Y2(I)/Y1(I)
           PSSUB(I)=PSDEP(I)
            PSDEP(I)=MAX(PSDEP(I), 0.)
            PSSUB(I)=MAX(-QS(I), MIN(PSSUB(I), 0.))
              Y2(I)=.78/ZG(I)**2+R20BQ*SCV(I)/ZG(I)**BGH5
           PGSUB(I)=R20T*SSI(I)*Y2(I)/Y1(I)
            DM(I)=QV(I)+QB0-QSI(I)
            RSUB1(I)=CS580*QSI(I)*RTAIR(I)*RTAIR(I)
C     ********   DEPOSITION OR SUBLIMATION OF QS  **********************
              Y1(I)=DM(I)/(1.+RSUB1(I))
           PSDEP(I)=MIN(PSDEP(I),MAX(Y1(I),0.))
              Y2(I)=MIN(Y1(I),0.)
           PSSUB(I)=MAX(PSSUB(I),Y2(I))
C     ********   SUBLIMATION OF QG   ***********************************
              DD(I)=MAX((-Y2(I)-QS(I)), 0.)
           PGSUB(I)=MIN(DD(I), QG(I), MAX(PGSUB(I),0.))
              DLT1(I)=CVMGP(1.,0.,QC(I)+QI(I)-1.E-5)
            PSDEP(I)=DLT1(I)*PSDEP(I)
            PSSUB(I)=(1.-DLT1(I))*PSSUB(I)
            PGSUB(I)=(1.-DLT1(I))*PGSUB(I)
           PT(I)=PT(I)+ASCP*(PSDEP(I)+PSSUB(I)-PGSUB(I))
           QV(I)=QV(I)+PGSUB(I)-PSSUB(I)-PSDEP(I)
           QS(I)=QS(I)+PSDEP(I)+PSSUB(I)
           QG(I)=QG(I)-PGSUB(I)
          ENDIF
C* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
         ERN(I)=0.0
          IF(QR(I).GT.0.0) THEN
            TAIR(I)=(PT(I)+TB0)*PI0
            RTAIR(I)=1./(TAIR(I)-C358)
            QSW(I)=RP0*EXP(C172-C409*RTAIR(I))
            SSW(I)=(QV(I)+QB0)/QSW(I)-1.0
            DM(I)=QV(I)+QB0-QSW(I)
            RSUB1(I)=CV409*QSW(I)*RTAIR(I)*RTAIR(I)
             DD1(I)=MAX(-DM(I)/(1.+RSUB1(I)), 0.0)
CM             Y1(I)=R00*QR(I)
CM         ERN(I)=(((1.6+124.9*Y1(I)**.2046)*Y1(I)**.525)
CM   1          /(2.55E6/(P00*QSW(I))+5.4E5))*(-DM(I)/(R00*QSW(I)))*D2T
             Y1(I)=.78/ZR(I)**2+R23AF*SCV(I)/ZR(I)**BWH5
             Y2(I)=R23BR/(TCA(I)*TAIR(I)**2)+1./(DWV(I)*QSW(I))
           ERN(I)=R23T*SSW(I)*Y1(I)/Y2(I)
           ERN(I)=MIN(DD1(I),QR(I),MAX(ERN(I),0.))
           PT(I)=PT(I)-AVCP*ERN(I)
           QV(I)=QV(I)+ERN(I)
           QR(I)=QR(I)-ERN(I)
          ENDIF
C        IF (QV(I)+QB0 .LE. 0.) QV(I)=-QB0
        IF (QC(I) .LE. cmin) QC(I)=0.
        IF (QR(I) .LE. cmin) QR(I)=0.
        IF (QI(I) .LE. cmin) QI(I)=0.
        IF (QS(I) .LE. cmin) QS(I)=0.
        IF (QG(I) .LE. cmin) QG(I)=0.
C        Y1(I)=PT(I)-XX0(I)
C        Y2(I)=QV(I)-XX00(I)
C       DPT(I,K)=PT(I)+TMD(I,K)
C       DQV(I,K)=QV(I)+QMD(I,K)
C        TMD(I,K)=Y1(I)
C        QMD(I,K)=Y2(I)
        DPT(I,K)=PT(I)
        DQV(I,K)=QV(I)
        QCL(I,K)=QC(I)
        QRN(I,K)=QR(I)
        QCI(I,K)=QI(I)
        QCS(I,K)=QS(I)
        QCG(I,K)=QG(I)
  275   CONTINUE
CC*********************************************************************
       SCC=0.
       SEE=0.
       DO 110 I=2,ILES
         DD(I)=MAX(-CND(I), 0.)
         CND(I)=MAX(CND(I), 0.)
         DD1(I)=MAX(-DEP(I), 0.)
         DEP(I)=MAX(DEP(I), 0.)
         SCC=SCC+CND(I)
  110    SEE=SEE+DD(I)+ERN(I)
        SC(K)=SCC+SC(K)
        SE(K)=SEE+SE(K)
CC    ***   STATISTICS FOR CONVECTIVE AND ANVIL REGIMES   ***********
       IF (ID .EQ. 1) THEN
c         RDTS=1./D2T
        DO 280 I=2,ILES
c           CND(I)=CND(I)*RDTS
c           DD(I)=DD(I)*RDTS
c           PINT(I)=PINT(I)*RDTS
c           PIDW(I)=PIDW(I)*RDTS
c           PIMLT(I)=PIMLT(I)*RDTS
c           PIHOM(I)=PIHOM(I)*RDTS
c           PSMLT(I)=PSMLT(I)*RDTS
c           PGMLT(I)=PGMLT(I)*RDTS
c           PSDEP(I)=PSDEP(I)*RDTS
c           PSSUB(I)=PSSUB(I)*RDTS
c           PGSUB(I)=PGSUB(I)*RDTS
c           DD1(I)=DD1(I)*RDTS
c           DEP(I)=DEP(I)*RDTS
c           ERN(I)=ERN(I)*RDTS
         DDA(I)=RSW(I,K)
         DDB(I)=RLW(I,K)
         Y1(I)=QC(I)+QR(I)+QI(I)+QS(I)+QG(I)
         dm(i)=a0*(rho1(k)*ww1(i,k)+rho1(k+1)*ww1(i,k+1)+
     1              Y0(I)*(rho1(k)*wb(k)+rho1(k+1)*wb(k+1)))
  280    RQ(I)=.005*(RHO1(K)*(WW1(I,K)+wb(k))+
     1               RHO1(K+1)*(WW1(I,K+1)+wb(k+1)))/r00


        DO 1050 KC=1,7
c        KC=4

         DO 320 MT=1,4
         DO 320 I=2,ILES
          IBZ(I,MT)=0
  320      IF(ICS(I,mt).EQ.1) IBZ(I,MT)=1
        DO 315 I=2,ILES
            IBZ(I,1)=1
  315    CONTINUE
         DO 330 MT=1,4
         DO 330 I=2,ILES
           IF(KC.EQ.4) GO TO 330
           IF(KC.LE.3) GO TO 36
            IF (RQ(I).GT.RBY(KC)) IBZ(I,MT)=0
            GO TO 330
   36       IF (RQ(I).LT.RBY(KC)) IBZ(I,MT)=0
  330     CONTINUE
         DO 350 MT=1,4
          SWW=0.0
          SCC=0.0
          SEE=0.0
          A1=0.0
          A2=0.0
          A3=0.0
          A4=0.0
          A5=0.0
          A6=0.0
          A7=0.0
          A8=0.0
          A9=0.0
          A10=0.0
          A11=0.0
          A12=0.0
          A13=0.0
          A14=0.0
          A15=0.0
          A16=0.0
          A17=0.0
          A18=0.0
          A19=0.0
          A20=0.0
          A21=0.0
          A22=0.0
          A23=0.0
          A24=0.0
          A25=0.0
          A26=0.0
          A27=0.0
          A28=0.0
          A29=0.0
          A30=0.0
          A31=0.0
          A32=0.0
          A33=0.0
          A34=0.0
          A35=0.0
          A36=0.0
          A37=0.0
          A38=0.0
          A39=0.0
          A40=0.0
          A41=0.0
          A42=0.0
          A43=0.0
          A44=0.0
          A45=0
          A46=0.
          A47=0.
          A48=0.
          A49=0.
          A50=0.

         DO 30 I=2,ILES
          IF(IBZ(I,MT).EQ.1) THEN
           SWW=SWW+DM(I)
           SCC=SCC+CND(I)
           SEE=SEE+DD(I)
           A1=A1+PINT(I)
           A2=A2+(PIDW(I)+PIHOM(I))
           A3=A3+PIMLT(I)
           A4=A4+PSAUT(I)
           A5=A5+PSACI(I)
           A6=A6+PSACW(I)
           A7=A7+PRACI(I)
           A8=A8+PIACR(I)
           A9=A9+PRAUT(I)
           A10=A10+PRACW(I)
           A11=A11+PSFW(I)
           A12=A12+PSFI(I)
           A13=A13+(PGACS(I)+DGACS(I)+WGACS(I))
           A14=A14+DGACW(I)
           A15=A15+(DGACI(I)+WGACI(I))
           A16=A16+DGACR(I)
           A17=A17+(WGACR(I)+WGACI(I)+WGACS(I)+DGACW(I))
           A18=A18+PGAUT(I)
           A19=A19+PRACS(I)
           A20=A20+PSACR(I)
           A21=A21+PGFR(I)
           A22=A22+PSMLT(I)
           A23=A23+PGMLT(I)
           A24=A24+PSDEP(I)
           A25=A25+PSSUB(I)
           A26=A26+PGSUB(I)
           A27=A27+DD1(I)
           A28=A28+PRACI(I)*DLT3(I)
           A29=A29+PIACR(I)*DLT3(I)
           A30=A30+PSACR(I)*DLT2(I)
           A31=A31+DEP(I)
           A32=A32+(DGACR(I)+DGACI(I)+DGACS(I)+DGACW(I))
           A33=A33+QC(I)
           A34=A34+QR(I)
           A35=A35+QI(I)
           A36=A36+QS(I)
           A37=A37+QG(I)
           A38=A38+ERN(I)
           A39=A39+WGACR(I)
           A40=A40+QSACW(I)
          A41=A41+DDA(I)
          A42=A42+DDB(I)
            A43=A43+(QV(I)+QA1(K)-QB(K))
            A44=A44+(PT(I)+TA1(K)-TB(K))
            A45=A45+1.
            A46=A46+Y1(I)
            A47=a47+(PSACW(I)+PSFW(I)+DGACW(I)+PIACR(I)+DGACR(I)
     1               +WGACR(I)+PSACR(I)+PGFR(I)-PIMLT(I)+PIDW(I))
           A48=A48+(Y1(I)-QCL1(I,K)-QRN1(I,K)-QCI1(I,K)-QCS1(I,K)
     1            -QCG1(I,K))
           A49=A49+(QV(I)-DQV1(I,K))
           A50=A50+(PT(I)-DPT1(I,K))

          ENDIF
   30   CONTINUE
        SMF0(K,MT,KC)=SWW+SMF0(K,MT,KC)
        COC(K,MT,KC)=SCC+COC(K,MT,KC)
        COE(K,MT,KC)=SEE+COE(K,MT,KC)
        THOM(K,MT,KC)=THOM(K,MT,KC)+A1
        TDW(K,MT,KC)=TDW(K,MT,KC)+A2
        TMLT(K,MT,KC)=TMLT(K,MT,KC)+A3
        SAUT(K,MT,KC)=SAUT(K,MT,KC)+A4
        SACI(K,MT,KC)=SACI(K,MT,KC)+A5
        SACW(K,MT,KC)=SACW(K,MT,KC)+A6
        RACI(K,MT,KC)=RACI(K,MT,KC)+A7
        TACR(K,MT,KC)=TACR(K,MT,KC)+A8
        RAUT(K,MT,KC)=RAUT(K,MT,KC)+A9
        RACW(K,MT,KC)=RACW(K,MT,KC)+A10
        SFW(K,MT,KC)=SFW(K,MT,KC)+A11
        SFI(K,MT,KC)=SFI(K,MT,KC)+A12
        GACS(K,MT,KC)=GACS(K,MT,KC)+A13
        GACW(K,MT,KC)=GACW(K,MT,KC)+A14
        GACI(K,MT,KC)=GACI(K,MT,KC)+A15
        GACR(K,MT,KC)=GACR(K,MT,KC)+A16
        GWET(K,MT,KC)=GWET(K,MT,KC)+A17
        GAUT(K,MT,KC)=GAUT(K,MT,KC)+A18
        RACS(K,MT,KC)=RACS(K,MT,KC)+A19
        SACR(K,MT,KC)=SACR(K,MT,KC)+A20
        GFR(K,MT,KC)=GFR(K,MT,KC)+A21
        SMLT(K,MT,KC)=SMLT(K,MT,KC)+A22
        GMLT(K,MT,KC)=GMLT(K,MT,KC)+A23
        SDEP(K,MT,KC)=SDEP(K,MT,KC)+A24
        SSUB(K,MT,KC)=SSUB(K,MT,KC)+A25
        GSUB(K,MT,KC)=GSUB(K,MT,KC)+A26
        PERN(K,MT,KC)=PERN(K,MT,KC)+A27
        D3RI(K,MT,KC)=D3RI(K,MT,KC)+A28
        D3IR(K,MT,KC)=D3IR(K,MT,KC)+A29
        D2SR(K,MT,KC)=D2SR(K,MT,KC)+A30
        D2RS(K,MT,KC)=D2RS(K,MT,KC)+A31
        GDRY(K,MT,KC)=GDRY(K,MT,KC)+A32
        ERNS(K,MT,KC)=ERNS(K,MT,KC)+A38
        WGRS(K,MT,KC)=WGRS(K,MT,KC)+A39
        QSWS(K,MT,KC)=QSWS(K,MT,KC)+A40
        SRSW(K,MT,KC)=SRSW(K,MT,KC)+A41
        SRLW(K,MT,KC)=SRLW(K,MT,KC)+A42

        QC0(K,MT,KC)=A33
        QR0(K,MT,KC)=A34
        QI0(K,MT,KC)=A35
        QS0(K,MT,KC)=A36
        QG0(K,MT,KC)=A37
        SQC0(K,MT,KC)=SQC0(K,MT,KC)+QC0(K,MT,KC)
        SQR0(K,MT,KC)=SQR0(K,MT,KC)+QR0(K,MT,KC)
        SQI0(K,MT,KC)=SQI0(K,MT,KC)+QI0(K,MT,KC)
        SQS0(K,MT,KC)=SQS0(K,MT,KC)+QS0(K,MT,KC)
        SQG0(K,MT,KC)=SQG0(K,MT,KC)+QG0(K,MT,KC)

        qv0l(k,mt,kc)=a43
        tt0l(k,mt,kc)=a44
        sgpt(k,mt,kc)=a45
        tsqq(k,mt,kc)=a46
        stt0l(k,mt,kc)=stt0l(k,mt,kc)+tt0l(k,mt,kc)
        sqv0l(k,mt,kc)=sqv0l(k,mt,kc)+qv0l(k,mt,kc)
        ssgpt(k,mt,kc)=ssgpt(k,mt,kc)+sgpt(k,mt,kc)
        tsqq1(k,mt,kc)=tsqq1(k,mt,kc)+tsqq(k,mt,kc)

        q1tl(k,mt,kc)=q1tl(k,mt,kc)+a47
        sqqdtl(k,mt,kc)=sqqdtl(k,mt,kc)+a48
        sqdtl(k,mt,kc)=sqdtl(k,mt,kc)+a49
        stdtl(k,mt,kc)=stdtl(k,mt,kc)+a50

  350   CONTINUE
 1050  CONTINUE
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
C    CONDENSATION:  cnd(ij)
C    EVAPORATION:   dd(ij)+ern(ij)
C    DEPOSITION:    dep(ij)+psdep(ij)+pint(ij)
C    SUBLIMATION:   dd1(ij)+pgsub(ij)-pssub(ij)
C    MELTING:       psmlt(ij)+pgmlt(ij)+pimlt(ij)
C    FREEZING:      pihom(ij)+pidw(ij)+psacw(ij)+psfw(ij)+dgacw(ij)+
c                   piacr(ij)+dgacr(ij)+wgacr(ij)+psacr(ij)+pgfr(ij)                
C    MASS FLUX:     DM(IJ)
C    CLOUD WATER:   QC(IJ)
C    RAIN:          
C    CLOUD ICE
C    SNOW
C    HAIL/GRAUPEL:
C----------------------------------------------------------------------
       DO 42 IJ=2,iles
        CND(IJ)=cnd(ij)*RFT
        ERN(IJ)=(ern(ij)+dd(ij))*RFT
        Y1(IJ)=(dep(ij)+psdep(ij)+pint(ij))*RFT
        Y2(IJ)=(dd1(ij)+pgsub(ij)-pssub(ij))*RFT
        Y3(IJ)=(psmlt(ij)+pgmlt(ij)+pimlt(ij))*RFT
        Y4(IJ)=(pihom(ij)+pidw(ij)+psacw(ij)+psfw(ij)+dgacw(ij)+
     1          piacr(ij)+dgacr(ij)+wgacr(ij)+psacr(ij)+pgfr(ij))*RFT
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
       DO 40 IJ=2,iles

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
 1000 CONTINUE

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
c        ELSE
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
      RETURN
      END
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TERVL (IRSG,RHO,FV)
C     COMPUTE THE TERMINAL VELOCITY OF QR, QS AND QG
      PARAMETER (NX=514,NZ=43)
      PARAMETER (NX10=10*NX)
      common/ilin/ lin

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC

      real    ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq

      real    qrn(nx,nz),qcg(nx,nz),qcs(nx,nz)
      common/b1r/ qrn
      common/b1g/ qcg
      common/b1s/ qcs

      real    ww1(nx,nz)
      COMMON/BTV/ WW1

      real    y1(nx),y2(nx),y3(nx),vr(nx),vs(nx),vg(nx),y4(nx10)
      common/ba/ y1,y2,y3,vr,vs,vg,y4

      integer irsg
      real    fv(nz),rho(nz)

C     ******************************************************************
      IF(IRSG.NE.0) GO TO 1
C     ***      IRSG=0    FOR RAIN WATER (QR)  **********
      DO 100 K=2,KLES
         VRCF=VRC*FV(K)
       DO 10 I=2,ILES
        WW1(I,K)=0.0
          A1=.5*RHO(K)*(QRN(I,K)+QRN(I,K-1))
        IF (A1 .GT. 1.E-15) THEN
         WW1(I,K)=VRCF*A1**BWQ
        ENDIF
   10  CONTINUE
  100 CONTINUE
      RETURN
CC    ***      IRSG=1    FOR SNOW (QS)        **********
    1 IF(IRSG.NE.1) GO TO 2
      DO 200 K=2,KLES
          VSCF=VSC*FV(K)
       DO 20 I=2,ILES
        WW1(I,K)=0.0
          A1=.5*RHO(K)*(QCS(I,K)+QCS(I,K-1))
         IF (A1 .GT. 1.E-15) THEN
          WW1(I,K)=VSCF*A1**BSQ
         ENDIF
   20  CONTINUE
  200 CONTINUE
      RETURN
CC    ***      IRSG=2    FOR GRAUPEL (QG)      **********
    2 DO 300 K=2,KLES
          IF (LIN .EQ. 1) THEN
            VGCR=VGC/SQRT(RHO(K))
          ELSE
            VGCR=VGC*FV(K)
          ENDIF
       DO 30 I=2,ILES
        WW1(I,K)=0.0
          A1=.5*RHO(K)*(QCG(I,K)+QCG(I,K-1))
         IF (A1 .GT. 1.E-15) THEN
          WW1(I,K)=VGCR*A1**BGQ
         ENDIF
   30  CONTINUE
  300 CONTINUE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sepca (isec,ics,iv,aco,aco1,aan,aan1,lconv,lanvl,lnspt)
cc   ***  gce model's convective and anvil separation   *****
      parameter (NX=514,NZ=43)
      parameter (nx13=16*nx-2*nx-nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    qc(nx,nz),qr(nx,nz),qi(nx,nz),qs(nx,nz),qg(nx,nz),w(nx,nz)

      common/b2c/ qc
      common/b2r/ qr
      common/b2i/ qi
      common/b2s/ qs
      common/b2g/ qg
      common/b2w/ w

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1

      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real y1(nx),y3(nx),tair(nz),y2(nx13)
      common/ba/ y1,y3,tair,y2

      integer ics(nx,4),iv(nx),jv(nx),it(nx),itt(nx)
cc
       do 4 i=1,imax
          ics(i,1)=1
    4  continue
c
       do 1 ist=2,4
       do 1 i=1,imax
    1    ics(i,ist)=0
       do 2 i=1,imax
         itt(i)=4
         if(ri(i).ge.0.001) itt(i)=1
         it(i)=4
         iv(i)=0
         jv(i)=0
    2    if (ri(i).lt.0.001) it(i)=1

      if(isec.eq.isec/14400*14400) then
        write(6,4567) 
 3456 format(/3x,128i1/)
 4567  format(/5x,'gce models convective-stratiform partition')
       write(6,3456) (it(i),i=2,iles,4)
       write(6,3456) (itt(i),i=2,iles,4)
      endif

      do 3 k=2,kles
    3   tair(k)=pi(k)*tb(k)-273.16
c   ***  churchill and houze's method (sfc rainfall only)  ********
       ri(1)=ri(iles)
       ri(imax)=ri(2)
        do 214 i=2,iles
         im=i-1
         imm=i-2
          if (i .eq. 2) imm=il2
         ip=i+1
         ipp=i+2
          if (i .eq. iles) ipp=3
  214   y1(i)=.4*(ri(imm)+ri(im)+ri(i)+ri(ip)+ri(ipp))
       nite=0
       do 245 i=2,iles
        if (ri(i) .ge. .001) nite=nite+1
  245  continue
        nstep=0
  248  nstep=nstep+1
        big=0.0
        ibig=2
       do 242 i=2,iles
        a1=ri(i)
        if (it(i) .le. 3) go to 242
        if (a1 .lt. big) go to 242
        ibig=i
        big=a1
  242  continue
        ip1=ibig+1
        im1=ibig-1
        rave=y1(ibig)
         if (ri(ibig) .lt. rave) go to 246
        k=3
        kp1=3
        km1=3
        go to 244
  246    k=2
         kp1=4
         km1=4
        if (it(ip1) .eq. 3) kp1=3
        if (it(im1) .eq. 3) km1=3
  244   it(ibig)=k
        if (it(ip1) .gt. 2) it(ip1)=kp1
        if (it(im1) .gt. 2) it(im1)=km1
       if (nstep .lt. nite) go to 248

      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif

c  ***  find the location and value of max sfc rainfall grid point
        big=0.0
        ibig=2
       do 251 i=2,iles
        a1=ri(i)
        if (a1 .lt. big) go to 251
        ibig=i
        big=a1
  251  continue

c        rianmax=max(20., min(25., big))
c        if(isec .le. 3*3600) rianmax=max(10., min(15.0, big))
c        if(isec .le. 6*3600) rianmax=max(15., min(20.0, big))
        rianmax=20.
       do 255 i=2,iles
         if (ri(i) .ge. rianmax) it(i)=3
  255    if (ri(i) .lt. 0.001) it(i)=1

      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif

c   ***  tao & simpson (1989) and tao et al (1993)   **************
c     ****   search for the largest value of w and cloud water
c      bigw=0.
c      do 10 k=2,kles
c      do 10 i=2,iles
c       xa=abs(w(i,k))
c       bigw=max (bigw,xa)
c   10 continue
c       wbig=max(300., min(300., 0.5*bigw), 0.25*bigw)
c      bigqc=0.
c      do 12 k=2,kles
c      do 12 i=2,iles
c       xa=qc(i,k)
c       bigqc=max (bigqc,xa)
c   12 continue
c       qcbig= max(0.50e-3, 0.5*bigqc)
c       qcbig1=max(1.00e-3, 0.5*bigqc)
       wbig=300.
       qcbig= 0.50e-3
       qcbig1=1.00e-3
cccccccccccccccccccccccccccccccccccccccccccccccc
         do 300 i=2,iles
          icloudw=0
          icloudl=0
          icloud=0
          icwwl=0
          do 30 k=2,kles
cc  ***   lower cloudy region   ****************
           if (tair(k) .ge. 0.0) then
            if (qc(i,k) .ge. qcbig) icloudw=1
            if (qc(i,k) .ge. qcbig1) icloudl=1
            if (w(i,k) .ge. wbig) icwwl=1
           else
cc  ***   middle and upper cloudy region
           endif
   30     continue
c   ***  middle-upper level w > 0.5*wbig m/s or low-level qc > 1 g/kg
          if (it(i).eq.1) then
            if (icwwl .eq. 1 .and. icloudw .eq. 1) it(i)=3
          endif
           if (it(i).eq.2) then
             if (icloudl .eq. 1) it(i)=3
           endif
  300   continue
c
         do 350 i=2,iles
          icloud=0
          do 35 k=2,kles
            rzm=rho(k)*1.e6
            y3(k)=rzm*(qc(i,k)+qi(i,k)+qr(i,k)+qs(i,k)+qg(i,k))
           if (y3(k) .ge. .01) icloud=1
   35     continue
          if (it(i).eq.1 .and. icloud.eq.1) jv(i)=1
  350   continue
cccccccccccccccccccccccccccccccccccccccccccccccc
      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif

c   ***   ********************************************************
c   ***   it(i)=3  convective region
c   ***   it(i)=2  stratiform region
c   ***   it(i)=4  stratiform region but no sfc precipitation
c   ***   it(i)=1  cloud free region

       a1=0.
       a2=0.
       a3=0.
       a4=0.
       lconv=0
       lanvl=0
       lnspt=0
       do 900 i=2,iles
        if (it(i) .eq. 2) then
          ics(i,3)=1
          a2=a2+ri(i)
          a4=a4+qr(i,2)
          lanvl=lanvl+1
        endif
        if (it(i) .eq. 3) then
          ics(i,2)=1
          a1=a1+ri(i)
          a3=a3+qr(i,2)
          lconv=lconv+1
        endif
        if (it(i) .eq. 1 .and. jv(i) .eq. 1) then
          ics(i,4)=1
          it(i)=4
          lnspt=lnspt+1
        endif
        iv(i)=it(i)
  900  continue
       aco=aco+a1
       aco1=aco1+a3
       aan=aan+a2
       aan1=aan1+a4
      if(isec.eq.isec/14400*14400) then
        write(6,3456) (it(i),i=2,iles,4)
      endif
      return
      end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tmaxadv (x,small)
C   ****   FIND MAXIMUM VALUES' routine
      PARAMETER (NX=514,NZ=43)
      real    X(NX,NZ)
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      save
      SMALL=0.
C     ****   SEARCH FOR THE SMALLEST VALUE
      DO 10 K=2,KLES
      DO 10 I=2,ILES
       XA=X(I,K)
       SMALL=MIN (SMALL,XA)
   10 CONTINUE
      RETURN
      END
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tmax (x,aa,aminut,big,small,id)
      implicit none
c   ****   find maximum values' routine
      integer nx,nz
      parameter (NX=514,NZ=43)
      character*4 aa
      integer id
      real    aminut,big,small
      real    x(nx,nz)
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    cons(4)
      data cons/.01,1.,1000.,.0000001/

      integer i,k,ilo1,ilo2,klo1,klo2
      real    xa


      save

c
      big=0.
      small=0.
c     ****   search for the largest or/and smallest values
      do 10 k=2,kles
      do 10 i=2,iles
       xa=x(i,k)
       big=max (big,xa)
       small=min (small,xa)
   10 continue
c     ****   search for the location
      do 20 k=2,kles
      do 20 i=2,iles
       if(x(i,k).ge.big) then
        ilo1=i-1
        klo1=k-1
       endif
       if(x(i,k).le.small) then
        ilo2=i-1
        klo2=k-1
       endif
   20 continue
c     ************************************
      big=cons(id)*big
      small=cons(id)*small
      write(6,123) aminut,aa,big,ilo1,klo1,small,ilo2,klo2
  123 format(3x,5hamin=,f6.0,3x,a4,4x,4hbig=,f7.3,3x,2hi=,i3,3x,2hk=,
     1   i3,6x,6hsmall=,f7.3,3x,2hi=,i3,2x,2hk=,i3)
      return
      end
C-----------------------------------------------------------------------------
      subroutine wtap (x,itape)
c     ******   write data to/from restart file   ******
      parameter (NX=514,NZ=43)
      real x(nx,nz)
      save
      if (itape.gt.6) go to 100
      write (itape) x
      return
  100  read (itape) x
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtapt (itape)
c     ******   write data to/from restart file  ******
      parameter (NX=514,NZ=43,NT=2880,ITT=244)
      parameter (nz2=2*nz,nt5=5*nt)
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     1   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     2   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/bb6/tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/b6/tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     $   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/brh1/z(nz2),qc(nz),qr(nz),qi(nz),qs(nz),qg(nz),sqc(nz),
     $        sqr(nz),sqi(nz),sqs(nz),sqg(nz),t(nt5)
      common/q1q2t/aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      COMMON/DINRAD/ P00(NZ),DZ0(NZ),TAIRSFC(NX),QAIRSFC(NX),THAIRSF(NX)
     *,       PAIRSFC(NX)
      common/sfcten/tsdt,qsdt,thsdt,psdt
      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     ****************************************
      if(itape.gt.6) go to 100
        write(itape) ub
        write(itape) vb
        write(itape) ub1
        write(itape) vb1
        write(itape) tb
        write(itape) qb
        write(itape) rho
        write(itape) rho1
        write(itape) ta
        write(itape) qa
        write(itape) ta1
        write(itape) qa1
        write(itape) pi
        write(itape) p0
        write(itape) wb
        write(itape) wbt
        write(itape) ubt
        write(itape) vbt ! 8/21/01
        write(itape) tls ! 8/21/01
        write(itape) qls ! 8/21/01
        write(itape) fd
        write(itape) fe
        write(itape) qc
        write(itape) qr
        write(itape) qi
        write(itape) qs
        write(itape) qg
        write(itape) sqc
        write(itape) sqr
        write(itape) sqi
        write(itape) sqs
        write(itape) sqg
        write(itape) q1t
        write(itape) q2t
        write(itape) aq1t
        write(itape) aq2t
        write(itape) aq1zt
        write(itape) aq2zt
        write(itape) tairsfc
        write(itape) qairsfc
        write(itape) thairsf
        write(itape) pairsfc
        write(itape) tsdt
        write(itape) qsdt
        write(itape) thsdt
        write(itape) psdt
      return
  100   read(itape) ub
        read(itape) vb
        read(itape) ub1
        read(itape) vb1
        read(itape) tb
        read(itape) qb
        read(itape) rho
        read(itape) rho1
        read(itape) ta
        read(itape) qa
        read(itape) ta1
        read(itape) qa1
        read(itape) pi
        read(itape) p0
        read(itape) wb
        read(itape) wbt
        read(itape) ubt
        read(itape) vbt ! 8/21/01
        read(itape) tls ! 8/21/01
        read(itape) qls ! 8/21/01
        read(itape) fd
        read(itape) fe
        read(itape) qc
        read(itape) qr
        read(itape) qi
        read(itape) qs
        read(itape) qg
        read(itape) sqc
        read(itape) sqr
        read(itape) sqi
        read(itape) sqs
        read(itape) sqg
        read(itape) q1t
        read(itape) q2t
        read(itape) aq1t
        read(itape) aq2t
        read(itape) aq1zt
        read(itape) aq2zt
        read(itape) tairsfc
        read(itape) qairsfc
        read(itape) thairsf
        read(itape) pairsfc
        read(itape) tsdt
        read(itape) qsdt
        read(itape) thsdt
        read(itape) psdt
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     subroutine wstrad (itape,cosz,month,iday,hrl)
      subroutine wstrad (itape,cosz,month,iday,hrl,rlat)  ! 8/21/01 shie
cccccccccccccccccc  write/read to restart file cccccccccccccccccccccc
      parameter (NX=514)
      real ri180(nx),riold180(nx)
      real cosz,hrl,rlat  ! 8/21/01 shie
      integer month,iday  ! 8/21/01 shie
      common/sfcri/ ri180,riold180,iricont
      save
      if(itape.gt.6) go to 100
       write(itape) month
       write(itape) iday
       write(itape) iricont
       write(itape) cosz
       write(itape) hrl
       write(itape) rlat  ! 8/21/01 shie
       write(itape) ri180
       write(itape) riold180
      return
  100  read(itape) month
       read(itape) iday
       read(itape) iricont
       read(itape) cosz
       read(itape) hrl
       read(itape) rlat  ! 8/21/01 shie
       read(itape) ri180
       read(itape) riold180
      return
      end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wstart (itape)
Ccccccccccccccc write to/from restart file cccccccccccccccccccccccccccccc
      parameter (NX=514,NZ=43,nnt=481)
      parameter (nz17=17*nz,nz4=4*nz,nz2=2*nz,nnt3=3*nnt,nx2=2*nx,
     1           nz9=9*nz,nnt7=7*nnt)
      common/bstart/ iw(3)
      common/rbstart/ rw(3)
c 8/21/01 shie found the dummy "sltm" was supposed to be as follows:
c  based on "gcem3d.F".  Put here for reference.
c COMMON/SLTM/ RLAT,RMONTH,RIDAY,HRL,SO0,COSZ,ICOSZ,TERMAN,RSFC
      common/sltm/ rlat(9) ! 8/21/01 shie, dummy!
      common/sun/ isun(2)  ! 8/21/01 shie, dummy!
      common/mbudget/ acoc(nz17)
      common/tqave/ tavet(nz2)
      common/tbudget/ avett(nnt3)
      common/tbudget1/ t_sfcq(nnt7)
      common/bsfc/ tsfc_1(nx2)
      common/gbs/ tlsw(nz9)
      common/gbs11/ tlsw1(nz4)
      common/timestat/ ndt_stat(3)
      COMMON/SFLUXS/ SUW(NX),SVW(NX),SWT(NX),SWQ(NX)
      save

      if(itape.gt.6) go to 100
        write(itape) iw
        write(itape) isun     ! 8/21/01 shie, dummy!
        write(itape) ndt_stat
        write(itape) rw
        write(itape) rlat     ! 8/21/01 shie, dummy!
c        write(itape) acoc
c        write(itape) tavet
        write(itape) avett
        write(itape) t_sfcq
c        write(itape) tsfc_1
c        write(itape) tlsw
c        write(itape) tlsw1
        write(itape) suw
        write(itape) svw
        write(itape) swt
        write(itape) swq

      return
  100   read(itape) iw
        read(itape) isun       ! 8/21/01 shie, dummy!
        read(itape) ndt_stat
        read(itape) rw
        read(itape) rlat       ! 8/21/01 shie, dummy!
c        read(itape) acoc
c        read(itape) tavet
        read(itape) avett
        read(itape) t_sfcq
c        read(itape) tsfc_1
c        read(itape) tlsw
c        read(itape) tlsw1
        read(itape) suw
        read(itape) svw
        read(itape) swt
        read(itape) swq
      return
      end

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wtap0(x,itape)
c     ******   write 1-d array data to data file   ******
      real x
      real*4 x4
      save
      x4=x
      write (itape) x4
      return
      end

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtap1(x,x4,n,itape)
c     ******   write 1-d array data to data file   ******
      real x(n)
      real*4 x4(n)
      save
      do k=1,n
        x4(k)=x(k)
      enddo
      write (itape) x4
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtap2(x,x4,n1,n2,itape)
c     ******   write 2-d array data to data file   ******
      real x(n1,n2)
      real*4 x4(n1,n2)
      save
      do k=1,n2
        do i=1,n1
          x4(i,k)=x(i,k)
        enddo
      enddo
      write (itape) x4
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wtap3(x,x4,n1,n2,n3,itape)
c     ******   write 3-d array data to data file   ******
      real x(n1,n2,n3)
      real*4 x4(n1,n2,n3)
      save
      do k=1,n3
        do j=1,n2
          do i=1,n1
            x4(i,j,k)=x(i,j,k)
          enddo
        enddo
      enddo
      write (itape) x4
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outrf  (ircl,npp1,pl,pa)
      implicit none
      integer lay,iopt
      parameter(lay=88,iopt=8)

      integer ircl,npp1
      real    pl(lay),pa(lay)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)
      common/srflx/sfir,sfsw,shir,shsw,salpha,si

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    fir(4),fsw(4),hsw(4),hir(4)
      integer k,lc

      save


      write (iopt,605) ircl,(si(lc),lc=1,4)
      write (iopt,610)
      do 200  k=1,npp1
      do 100 lc=1,4
      if(si(lc).eq.0.0) go to 100
          fir(lc)=sfir(k,lc)/si(lc)
          fsw(lc)=sfsw(k,lc)/si(lc)
          hir(lc)=shir(k,lc)/si(lc)
          hsw(lc)=shsw(k,lc)/si(lc)
  100 continue
        write (iopt,600) k,pl(k),(fir(lc),lc=1,4),(fsw(lc),lc=1,4),
     +                     pa(k),(hir(lc),lc=1,4),(hsw(lc),lc=1,4)
  200 continue
      do 300 lc=1,4
      salpha(lc)=salpha(lc)/ircl
  300 continue
      write (iopt,615) salpha
  600 format(1x,i2,f8.2,8f8.2,f8.2,8f6.2)
  605 format(1x,'time count=',i4,10x,'cloud count=',4f10.0)
  610 format(2x,'k',6x,'pl',2x,'fi-tot',3x,'fi-cv',3x,'fi-cm',2x,
     +'fi-clr',2x,'fs-tot',3x,'fs-cv',3x,'fs-cm',2x,'fs-clr',
     + 4x,'pa',2x,'hi-tot',1x,'hi-cv',1x,'hi-cm',   'hi-clr',
     +            'hs-tot',1x,'hs-cv',1x,'hs-cm',   'hs-clr')
  615 format(1x,'albedo=',4f8.5)

      return
      end
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
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine humidity(T,P,Qsat)                                 
c
c     Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 
c
      Qsat = (1.0007+3.46e-6*P)*6.1121*exp(17.502*T/(240.97+T))
      return
      end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LKB(RR,RT,IFLAG)
C
C       TO DETERMINE THE LOWER BOUNDARY VALUE RT OF THE 
C       LOGARITHMIC PROFILES OF TEMPERATURE (IFLAG=1) 
C       OR HUMIDITY (IFLAG=2) IN THE ATMOSPHERE FROM ROUGHNESS 
C       REYNOLD NUMBER RR BETWEEN 0 AND 1000.  OUT OF RANGE
C       RR INDICATED BY RT=-999. BASED ON LIU ET AL.(1979)
C       JAS 36 1722-1723
c     New scalar RR relation from Moana Wave data.
C
      DIMENSION A(9,2),B(9,2),RAN(9)
      DATA A/0.177,2.7e3,1.03,1.026,1.625,4.661,34.904,1667.19,5.88E5,
     10.292,3.7e3,1.4,1.393,1.956,4.994,30.709,1448.68,2.98E5/
      DATA B/0.,4.28,0,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,
     10.,4.28,0,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/
      DATA RAN/0.11,.16,1.00,3.0,10.0,30.0,100.,300.,1000./
      I=1
      IF(RR.LE.0..OR.RR.GE.1000.) GOTO 90
   10 CONTINUE
      IF(RR.LE.RAN(I)) GOTO 20
      I=I+1
      GOTO 10
   20 RT=A(I,IFLAG)*RR**B(I,IFLAG)
      GOTO 99
   90 RT=-999.
   99 RETURN
      END
c
c
      FUNCTION PSI(ID,ZL)
C
C       TO EVALUATE THE STABILITY FUNCTION PSI FOR WIND SPEED (IFLAG=1)
C       OR FOR TEMPERATURE AND HUMIDITY PROFILES FROM STABILITY 
C       PARAMETER ZL. SEE LIU ET AL (1979).
c       Modified to include convective form following Fairall (Unpublished)
C
      IF(ZL)10,20,30                                                 
   10   F=1./(1+zl*zl) 
        CHIK=(1.-16.*ZL)**0.25
         IF(ID.EQ.1) GOTO 11
        PSIK=2.*LOG((1.+CHIK*CHIK)/2.)
        GOTO 12
   11   PSIK=2.*LOG((1.+CHIK)/2.)+LOG((1.+CHIK*CHIK)/2.)
     1 -2.*ATAN(CHIK)+2.*ATAN(1.)
   12   CHIC=(1.-12.87*ZL)**.333    !for very unstable conditions
        PSIC=1.5*LOG((CHIC*CHIC+CHIC+1.)/3.)
     &      -(3.**.5)*ATAN((2*CHIC+1.)/(3.**.5))
     &      +4.*ATAN(1.)/(3.**0.5) 
c                     
c     match Kansas and free-conv. forms with weighting F
c
      PSI= F*PSIK+(1-F)*PSIC
        goto 99
   20 PSI=0.
      GOTO 99
   30 continue
       PSI=-4.7*ZL
   99 RETURN
      END
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ZETA(T,Q,USR,TSR,QSR,Z,ZL)
C
C       TO EVALUATE OBUKHOVS STABILITY PARAMETER Z/L FROM AVERAGE
C       TEMP T IN DEG C, AVERAGE HUMIDITY Q IN GM/GM, HEIGHT IN M,
C       AND FRICTIONAL VEL,TEMP.,HUM. IN MKS UNITS
C       SEE LIU ET AL. (1979)
C
      common/const/al,beta,cpa,cpw,grav,hl,rhoa,rhow,rgas,toK,
     &visa,visw,von,fdg
c
      TA=T+toK
      TV=TA*(1.+0.61*Q)
      TVSR=TSR*(1.+0.61*Q)+0.61*TA*QSR
      IF(TVSR.EQ.0.)GOTO 10
      OB=TV*USR*USR/(grav*VON*TVSR)
      ZL=Z/OB
      GOTO 99
   10 ZL=0.
   99 RETURN
      END
       
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

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ptpb (m1)
      parameter (NX=514,NZ=43,NT=2880,nxi=nx-2,nzm1=nz-1)
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1 st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b1c/ qcl(nx,nz)
      common/b1r/ qrn(nx,nz)
      common/b1i/ qci(nx,nz)
      common/b1s/ qcs(nx,nz)
      common/b1g/ qcg(nx,nz)
      common/bw1/ rst(nt,nxi),pcltop(nt,nxi),pclbot(nt,nxi)
      dimension tqe(nz)
      save
      do 100 i=2,nx-1
      ix=i-1
       do k=2,nz-1
         tqe(k)=qcl(i,k)+qrn(i,k)+qci(i,k)+qcs(i,k)+qcg(i,k)
       end do
         call ttop(tqe,ktop,kbase)
          pcltop(m1,ix)=1.e-3*p0(ktop+1)
            if ( ktop.eq.nzm1) pcltop(m1,ix)=1.e-3*p0(ktop)
          pclbot(m1,ix)=1.e-3*p0(kbase)
            if (kbase.eq. 2  ) pclbot(m1,ix)=1.e-3*p0(kbase)
           thin=pclbot(m1,ix)-pcltop(m1,ix)
            if(thin.lt.50.) go to 100
  100 continue
      return
      end          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ttop (twc,ktop,kbase)

      implicit none
      integer nz,nzp,nzm
      parameter (NZ=43,nzp=nz+1,nzm=nz-1)

      integer kbase,ktop
      real    twc(nz)

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer k,k1,ktemp

      save


c
        ktop=2
        do 10 k=2,nzm
         if (twc(k).lt.1.e-5) goto 10
          ktemp=k
         if(ktop.lt.ktemp) ktop=ktemp
   10   continue
         kbase=nzm
        do 20 k1=2,nzm
         k=nzp-k1
         if (twc(k).lt.1.e-5) goto 20
          ktemp=k
          if (kbase.gt.ktemp) kbase=ktemp
   20    continue
      if(ktop.lt.kbase)then
        ktop=1
        kbase=2
      endif
      return
      end   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zang (month,iday,cosz,rlat,hrl)

      implicit none
      integer month,iday
      real    cosz,rlat,hrl

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer imd(12)
      data imd/31,28,31,30,31,30,31,31,30,31,30,31/

      integer idayi,monthi
      real    ang,anx,cosphi,costph,costtp,delta,eq,fac
      real    phi,sinphi,sintph,sinttp,tphi,ttphi,x,xlat

      save

c
c      noczon=0
c      fixmon=12.
      if (month .gt. 12) month=month-12 
      if (hrl .ge. 24.) then
       hrl=hrl-24.0
       iday=iday+1
         if (iday .gt. imd(month)) then
          iday=iday-imd(month)
          month=month+1
         end if
      end if
       monthi=month
        idayi=iday
c      if (noczon .eq. 1) then
c      monthi=fixmon     
c      end if
       fac=.01745
       x=float((monthi-1)*30+idayi+monthi/2)
       xlat=rlat*fac
       print*,'month day   hour     cosz'
       phi=2.*3.14159*x/365.
       tphi=phi*2.
       ttphi=phi*3.
       cosphi=cos(phi)
       sinphi=sin(phi)
       costph=cos(tphi)
       sintph=sin(tphi)
       costtp=cos(ttphi)
       sinttp=sin(ttphi)
       delta=.006918-0.399912*cosphi+0.070257*sinphi-0.006758*
     *        costph+0.000907*sintph-0.002697*costtp+.00148*sinttp
       anx=(hrl-12.)*15.*fac
       eq=.000075+.001868*cosphi-.032077*sinphi-.014615*costph
     *     -.040849*sintph
       ang=anx+eq
       cosz=sin(xlat)*sin(delta)+cos(xlat)*cos(delta)*cos(ang)
      write (6,200) monthi,idayi,hrl,cosz
  200 format(1x,2i4,2f10.4)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine radrat (iflag,cosz,npp1,pl,pa)
c nb is a number of blocks for calculation
      integer nbb,nq,nx,nz,lay,nadd,nxi,nw,nx1,nx2,nz2,nz4,nz15,nx3,mb1
      integer mb2,iflag,npp1
      parameter(NX=514,NZ=43,lay=88,nadd=7,nbb=1)
      parameter(nq=nz+nadd-1,nxi=nx-2,nw=nq-1,nx1=nxI/nbb/2,nx2=nxi/nbb)
      parameter(nz2=nz*2,nz4=nz*4,nz15=nz*15,nx3=nx*3)
      parameter(mb1=nz15,mb2=nz*7+nx3)

      real    cosz,pl(lay),pa(lay)
      integer iradave
      common/radflux/rflux(nxi/2,7)
      common/iptionr/ iradave
      common/iceopt/ ice913,ilif
      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      real    rsw(nx,nz),rlw(nx,nz)
      common/slwave/ rsw,rlw
      real    tb(nz),qb(nz),zz0(nz),rho(nz),zz1(nz2),tb1(nz),
     $        qb1(nz),zz2(mb1)
      common/b5/ tb,qb,zz0,rho,zz1,tb1,qb1,zz2
      real    zz3(nz4),p0(nz),pi(nz),zz4(mb2)
      common/b6/ zz3,p0,pi,zz4
      real    p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),pairsfc(nx),
     1        thairsf(nx)
      common/dinrad/ p00,dz0,tairsfc,qairsfc,thairsf,pairsfc
      real    rsirbm(nx1,1),rsirdf(nx1,1),rsuvbm(nx1,1),
     $        rsuvdf(nx1,1),taual(nx1,1,nw)
      common/albedo1/ rsirbm,rsirdf,rsuvbm,rsuvdf,taual
      integer ict,icb
      common/cloudp/ ict,icb
      real    sun_4(nx,4)
      common/surface/ sun_4
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    ta(lay),wa(lay),oa(lay),plu(nadd)
      real    fcld(nx1,1,nw),tauir1(nx1,1,nw),flx(nx1,1,nq)
      real    flc(nx1,1,nq),dfdts(nx1,1,nq),st4(nx1,1),tsq(nx1,1)
      real    plq(nx1,1,nq),taq(nx1,1,nw),waq(nx1,1,nw),oaq(nx1,1,nw)
      real    coolr1(nx,1,nw),heatr1(nx,1,nw)
      real    tausw1(nx1,1,nw,2),reff(nx1,1,nw,2)
      real    flx1(nx1,1,nq)
      real    cosz1(nx1,1)
      real    fdirir(nx1,1),fdifir(nx1,1)
      real    fdirpar(nx1,1),fdifpar(nx1,1)
      logical high
      data high/.false./
      data plu/ 0.01, 0.56, 3.04, 4.04, 7.18,  12.38, 22.46/
      integer i,k,irst,np,km,ib,i1,ii,i3,k1,ihalfp
      real    aco,cmc,sc0
      DATA IRST/0/
      save

      do k=1,nq
      do i=1,nx1
        flx(i,1,k)=0.0
        flx1(i,1,k)=0.0
        flc(i,1,k)=0.0
        dfdts(i,1,k)=0.0
      enddo
      enddo
      do k=1,nw
      do i=1,nx
        coolr1(i,1,k)=0.0
        heatr1(i,1,k)=0.0
      enddo
      enddo
      do i=1,nx1
        st4(i,1)=0.0
        fdirir(i,1)=0.0
        fdifir(i,1)=0.0
        fdirpar(i,1)=0.0
        fdifpar(i,1)=0.0
      enddo

      do i=1,lay
        pa(i)=0.0
      enddo
c
      IF (IRST.NE.0) GO TO 500
c
c     np=number of atmospheric layers; npp1=surface level
      np  =kl2+nadd
      npp1=np+1
c      cpi =4.*atan(1.)
c  solar constant and cosine of solar zenith angle
c      sc0=1365.
      sc0=1411.
      if (ilif .eq. 1) sc0=1365.
c     cosz=cos(51.74*cpi/180.)
        do i=1,nx1
          rsirbm(i,1)=0.07
          rsirdf(i,1)=0.07
          rsuvbm(i,1)=0.07
          rsuvdf(i,1)=0.07
        enddo
c
        ict=24
        icb=32
c 
c  assign co2 (cmc). units are parts/part
      cmc=300.e-6
      write(6,*) 'npp1(3)=',npp1
c     write(6,*) 'npp1(3)=',npp1,'  nadd = ',nadd
c  assign aerosol optical thickness 
      do 5 k=1,nw
      do 5 i=1,nx1
        taual(i,1,k)=0.0
    5 continue

      do k=1,kl2+1
        km=kl2+3-k
        pl(k+nadd)=1.e-3*p00(km)
      end do

      do k=1,nadd
        pl(k)=plu(k)
      end do

      do k=1,kl2
        km=kl2+2-k
        pa(k+nadd)=1.e-3*p0(km)
      end do

      do 13 k=1,nadd
        pa(k)=0.5*(pl(k)+pl(k+1))
        wa(k)=1.e-6
c        ta(k)=190.
   13 continue
c     print *,'       radrat.f          5  fito3.f'
      call fito3 (npp1,pa,ta,oa)
       print*
      print*,' k       pa        pl       taual        ao        ta
     *   wa'
      write (6,975) (k,pa(k),pl(k),taual(1,1,k),oa(k)
     1 ,ta(k),wa(k),k=1,np)    ! shie 2/22/01
c    1 ,ta(k),wa(k),k=1,npp1)

      k=npp1                   ! shie 2/22/01
      write (6,975) k,pa(k),pl(k),taual(1,1,np),oa(k) ! shie 2/22/01
     1 ,ta(k),wa(k)    ! shie 2/22/01
c
      IRST=1
c
  500 continue
  975 format(i3,3f10.3,e12.4,f12.5,f12.7)
        do i=1,nx1
          cosz1(i,1)=cosz
        enddo

      DO 1000 IB=1,NBB
        I1=(IB-1)*NX2+1
        CALL OPT4 (I1,NX1,NX2,PL,TA,WA,OA,TAUSW1,TAUIR1
     +                 ,FCLD,TAQ,WAQ,OAQ,PLQ,TSQ,REFF)

C ----------------- GCSS -----------------------------------

        do i=1,nx1
          rflux(i,1)=0.
          rflux(i,2)=0.
          rflux(i,3)=0.
          rflux(i,4)=0.
          rflux(i,5)=0.
          rflux(i,6)=0.
          rflux(i,7)=0.
        enddo

C ------------------------------------------------------------

        if (cosz.ge.0.005) then

          call sorad (nx1,1,1,plq,taq,waq,oaq,cmc,
     $                 tausw1,reff,taual,rsirbm,rsuvbm,
     1                 cosz1,flx1,fdirir,fdifir,fdirpar,fdifpar)
        end if

C ----------------- GCSS -----------------------------------

        do i=1,nx1
          rflux(i,2)=0.
        enddo

C ------------------------------------------------------------

        call irrad (nx1,1,1,tauir1,fcld,plq,taq,waq,oaq,cmc,tsq,
     $                  high,flx,flc,dfdts,st4)
c
c since flx1 is normalized, they should be multiplied
c by sc0*cosz
c since upward flux should be positive, heatr is multiplied by
c a minus sign
c save solar and long-wave radiative fluxes
C
C
        DO I=I1,I1+NX2-1,2
          II=1
          IF (IRADAVE .EQ. 0) II=(I-I1)/2+1
           i3=i+1
          sun_4(i3,1)=flx1(ii,1,nw+1)*sc0*cosz1(ii,1)
          sun_4(i3,2)=flx(ii,1,nw+1)
        enddo

C
        do k=1,nw
        DO I=I1,I1+NX2-1,2
          II=1
          IF (IRADAVE .EQ. 0) II=(I-I1)/2+1
           I3=I+1
          heatr1(i3,1,k)=(flx1(ii,1,k+1)-flx1(ii,1,k))*8.441874/
     1       (plq(ii,1,k+1)-plq(ii,1,k))
          heatr1(i3,1,k)=-heatr1(i3,1,k)*sc0*cosz1(ii,1)
          coolr1(i3,1,k)=(flx(ii,1,k+1)-flx(ii,1,k))*8.441874/
     1      (plq(ii,1,k+1)-plq(ii,1,k))
        enddo
        enddo
C
        do k1=nadd+1,np
          k=np+2-k1
        DO I=I1,I1+NX2-1,2
            i3=i+1
           aco=1./(86400.*pi(k))
          rsw(i3,k)=aco*heatr1(i3,1,k1)
          rlw(i3,k)=aco*coolr1(i3,1,k1)
        enddo
        enddo

1000   continue
CC
        IF (IRADAVE .EQ. 0) THEN
C
        DO I=2,NX-2,2
            i3=i+1
          sun_4(i3,1)=0.5*(sun_4(i3-1,1)+sun_4(i3+1,1))
          sun_4(i3,2)=0.5*(sun_4(i3-1,2)+sun_4(i3+1,2))
        enddo
        sun_4(nx-1,1)=.5*(sun_4(nx-2,1)+sun_4(2,1))
        sun_4(nx-1,2)=.5*(sun_4(nx-1,2)+sun_4(2,2))
        do k1=nadd+1,np
          k=np+2-k1
        DO I=2,NX-2,2
           i3=i+1
          rsw(i3,k)=0.5*(rsw(i3-1,k)+rsw(i3+1,k))
          rlw(i3,k)=0.5*(rlw(i3-1,k)+rlw(i3+1,k))
        enddo
          rsw(nx-1,k)=.5*(rsw(nx-2,k)+rsw(2,k))
          rlw(nx-1,k)=.5*(rlw(nx-2,k)+rlw(2,k))
        enddo
CC
       ELSE
C
        DO I=2,NX-2,2
            I3=I+1
          SUN_4(I3,1)=SUN_4(2,1)
          SUN_4(I3,2)=SUN_4(2,2)
        ENDDO
        SUN_4(NX-1,1)=SUN_4(2,1)
        SUN_4(NX-1,2)=SUN_4(2,2)
        DO K1=NADD+1,NP
          K=NP+2-K1
        DO I=2,NX-2,2
           I3=I+1
          RSW(I3,K)=RSW(2,K)
          RLW(I3,K)=RLW(2,K)
        ENDDO
          RSW(NX-1,K)=RSW(2,K)
          RLW(NX-1,K)=RLW(2,K)
        ENDDO
       ENDIF

c ------------------------- for GCSS workshop ----------------------------

       do i=1,nx1
         rflux(i,1)=rflux(i,1)*sc0*cosz1(i,1)
         rflux(i,3)=rflux(i,3)*sc0*cosz1(i,1)
c         rflux(i,3)=1411*sc0*cosz1(i,1)
         rflux(i,4)=flx1(i,1,1)*sc0*cosz1(i,1)-rflux(i,3)
         rflux(i,6)=flx1(i,1,npp1)*sc0*cosz1(i,1)-rflux(i,1)
       enddo


c -----------------------------------------------------------------------

      if(iflag.eq.2) then
        ihalfp=il2/2+1
        write(6,*) 'cosz=',cosz
        write(6,*) 'check rsw(i,k) and rlw(i,k) at central point'
        DO K1=NADD+1,NP
          K=NP+2-K1
          write(6,7821) k,pl(k),rsw(ihalfp,k),rlw(ihalfp,k),
     1        heatr1(ihalfp,1,k1),coolr1(ihalfp,1,k1)
        enddo
      endif
7821  format(2x,i6,2x,f10.3,2x,4e20.10)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine opt4(i1,nx1,nx2,pl,ta,wa,oa,tausw,tauir
     +               ,fcld,taq,waq,oaq,plq,tsq,reff)
c define variables and calculate the optical thickness
      integer nx,nz,lay,nadd,nz2,nz4,nz15,nx3
      integer npp1,mb1,mb2,nq,nw
      parameter(NX=514,NZ=43,lay=88)
      parameter(nadd=7,nz2=nz*2,nz4=nz*4,nz15=nz*15,nx3=nx*3)
      parameter(npp1=nz+nadd-1,mb1=nz15,mb2=nz*7+nx3,nq=npp1,nw=nq-1)

      integer i1,nx1,nx2
      real    pl(lay),ta(lay),wa(lay),oa(lay)
      real    tausw(nx1,1,nw,2),tauir(nx1,1,nw),reff(nx1,1,nw,2)
      real    fcld(nx1,1,nw),tsq(nx1,1)
      real    plq(nx1,1,nq),taq(nx1,1,nw),waq(nx1,1,nw),oaq(nx1,1,nw)

      integer iradave
      common/iptionr/ iradave
      COMMON/IPTIONR1/ IOPCLOUD
      common/iceopt/ ice913,ilif
      real    tnw,tns,tng,roqs,roqg,roqr
      common/size/ tnw,tns,tng,roqs,roqg,roqr

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),
     1        thairsf(nx),pairsfc(nx)
      common/dinrad/ p00,dz0,tairsfc,qairsfc,thairsf,pairsfc
      real    dpt(nx,nz),dqv(nx,nz),qci(nx,nz),qcl(nx,nz),
     $        qrn(nx,nz),qcs(nx,nz),qcg(nx,nz)
      common/b1t/ dpt
      common/b1q/ dqv
      common/b1i/ qci
      common/b1c/ qcl
      common/b1r/ qrn
      common/b1s/ qcs
      common/b1g/ qcg
      real    tb(nz),qb(nz),zz0(nz),rho(nz),zz1(nz2),tb1(nz),
     $        qb1(nz),zz2(mb1)
      common/b5/ tb,qb,zz0,rho,zz1,tb1,qb1,zz2
      real    zz3(nz4),p0(nz),pi(nz),zz4(mb2)
      common/b6/ zz3,p0,pi,zz4
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
c     real    qci(nx1,nz),qcl(nx1,nz),qrn(nx1,nz),qcs(nx1,nz)
c     real    aa1(nx,npp1),bb(nx,npp1),qcg(nx1,nz)
      real    aa1(nx,npp1),bb(nx,npp1),aa2(nx,npp1)
      real    y1(nw),y2(nw),y3(nw),y4(nw),y5(nw),y6(nw),y7(nw)
      integer i,i2,i3,ii,k,km
      real    b0,b1,b2,b4,b5,b6,b7,b9,cpi,effrad
      real    q1,q2,q3,q5,rnx1,tauqc,tauqg,tauqii,tauqis,tauqr,twco
      save
c
c      twcz=1.e-5
      I2=I1+NX2-1
      DO I=I1,I2,2
        II=(I-I1)/2+1
        I3=I+1
        do k=1,npp1
          plq(ii,1,k)=pl(k)
        enddo
        do k=1,nw
          oaq(ii,1,k)=oa(k)
          fcld(ii,1,k)=1.0
        enddo
        do k=1,nadd
          waq(ii,1,k)=wa(k)
          taq(ii,1,k)=ta(k)
c          taq(ii,1,k)=190.0
        enddo
        do k=2,kles
           km=kmax-k+nadd
          WAQ(II,1,KM)=DQV(I3,K)+QB1(K)
            IF(WAQ(II,1,KM).LE.1.E-6) WAQ(II,1,KM)=1.E-6
          TAQ(II,1,KM)=(DPT(I3,K)+TB1(K))*PI(K)
            IF(TAQ(II,1,KM).LE.190.) TAQ(II,1,KM)=190.
        enddo
        tsq(ii,1)=tairsfc(i3)
         if (ilif .eq. 1) tsq(ii,1)=29.+273.16
      enddo
       cpi =4.*atan(1.)
c       twcz=1.e-5
       twco=1.e-6
       DO 200 I=1,NX1
          II=I1+(I-1)*2
          i3=ii+1
         do k=1,nw
           tausw(i,1,k,1)=0.
           tausw(i,1,k,2)=0.
           reff(i,1,k,1)=0.
           reff(i,1,k,2)=0.
           tauir(i,1,k)=0.
         enddo
c        do k=2,kles
c        if (tqe(i,k) .ge. twcz) go to 30
c        enddo
c        go to 200
c30      continue
         b5=0.
         b6=0.
         b7=0.
         b9=0.
c         ttauqc=0.
c         ttauqr=0.
c         ttauis=0.
c         ttauii=0.
c         ttauqg=0.
         do 100 k=2,kles
           km=kmax-k+nadd
c         tausw(i,1,km)=0.
c         tauir(i,1,km)=0.
            b0=0.
            b1=0.
            b2=0.
            b4=0.
            b3=0.
            tauqc=0.0
            tauqr=0.0
            tauqs=0.0
            tauqis=0.0
            tauqii=0.0
            tauqg=0.0
           q1=qcl(i3,k)
           q2=qrn(i3,k)
           IF (IOPCLOUD .EQ. 1) THEN
              Q3=QCI(I3,K)+QCS(I3,K)
              Q4=0.
           ELSE
              Q3=QCI(I3,K)
              Q4=QCS(I3,K)
           ENDIF
           q5=qcg(i3,k)
          if(q1 .ge. twco) then
            b0=rho(k)*dz0(k)*q1
           effrad=0.0015
           tauqc=b0/effrad
            b0=1.e4*b0
            reff(i,1,km,2)=effrad*1.0e4
          endif
          if(q2 .ge. twco) then
            b1=rho(k)*dz0(k)*q2
           effrad=3./((cpi*tnw*roqr/(rho(k)*q2))**.25)
           tauqr=b1/effrad
            b1=1.e4*b1
          endif
          if(q4 .ge. twco) then
            b3=rho(k)*dz0(k)*q4
           effrad=3./((cpi*tns*roqs/(rho(k)*q4))**.25)
           tauqs=b3/effrad
            b3=1.e4*b3
          endif
          if(q5 .ge. twco) then
            b4=rho(k)*dz0(k)*q5
           effrad=3./((cpi*tng*roqg/(rho(k)*q5))**.25)
           tauqg=b4/effrad
            b4=1.e4*b4
          endif
            b4=b4+b3
            tauqg=tauqg+tauqs
          if(q3 .ge. twco) then
            b2=1.e4*rho(k)*dz0(k)*q3
c          effrad=0.0050
           effrad=0.0125+(taq(i,1,km)-243.16)*0.00050
           if (taq(i,1,km) .gt. 243.16) effrad=0.0125
           if (taq(i,1,km) .lt. 223.16) effrad=0.0025
           tauqis=b2*(-0.006656+ 3.686e-4/effrad)
           tauqii=b2*(-0.011500+ 4.110e-4/effrad
     +                         +17.300e-8/(effrad*effrad))
          reff(i,1,km,1)=effrad*1.0e4
          endif
          b5=b5+b0
          b6=b6+b1
          b7=b7+b2
          b9=b9+b4
c           ttauqc=ttauqc+tauqc
c           ttauqr=ttauqr+tauqr
c           ttauis=ttauis+tauqis
c           ttauii=ttauii+tauqii
c           ttauqg=ttauqg+tauqg
c          tausw(i,1,km)=1.5*(tauqc+tauqr)
           tausw(i,1,km,2)=1.5*(tauqc+tauqr)
           tauir(i,1,km)=0.5*tausw(i,1,km,2)
c          tauir(i,1,km)=0.5*tausw(i,1,km)
c          tausw(i,1,km)=tausw(i,1,km)+tauqis+tauqg
           tausw(i,1,km,1)=tauqis+tauqg
           tauir(i,1,km)=tauir(i,1,km)+tauqii+tauqg
  100    continue
  200    continue
         DO I=1,NX1
           k=2
           km=kmax-k+nadd
           k=3
           km1=kmax-k+nadd
           aa1(i,km)=0.8*tausw(i,1,km,1)+0.2*tausw(i,1,km1,1)
           aa2(i,km)=0.8*tausw(i,1,km,2)+0.2*tausw(i,1,km1,2)
           bb(i,km)=0.8*tauir(i,1,km)+0.2*tauir(i,1,km1)
cc           aa1(i,km)=tausw(i,1,km,1)
cc           aa2(i,km)=tausw(i,1,km,2)
cc           bb(i,km)=tauir(i,1,km)
           do 110 k=3,kles-1
             km=kmax-k+nadd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            aa1(i,km)=0.2*tausw(i,1,km-1,1)+0.60*tausw(i,1,km,1)
     1                +0.2*tausw(i,1,km+1,1)
             aa2(i,km)=0.2*tausw(i,1,km-1,2)+0.60*tausw(i,1,km,2)
     1                +0.2*tausw(i,1,km+1,2)
             bb(i,km)=0.2*tauir(i,1,km-1)+0.60*tauir(i,1,km)
     1               +0.2*tauir(i,1,km+1)
cc             aa1(i,km)=tausw(i,1,km,1)
cc             aa2(i,km)=tausw(i,1,km,2)
cc             bb(i,km)=tauir(i,1,km)
  110      continue
           k=kles
           km=kmax-k+nadd
           k=kl2
           km1=kmax-k+nadd
            aa1(i,km)=0.8*tausw(i,1,km,1)+0.2*tausw(i,1,km1,1)
            aa2(i,km)=0.8*tausw(i,1,km,2)+0.2*tausw(i,1,km1,2)
            bb(i,km)=0.8*tauir(i,1,km)+0.2*tauir(i,1,km1)
cc           aa1(i,km)=tausw(i,1,km,1)
cc           aa2(i,km)=tausw(i,1,km,2)
cc           bb(i,km)=tauir(i,1,km)
           do 120 k=2,kles
              km=kmax-k+nadd
             tausw(i,1,km,1)=aa1(i,km)
             tausw(i,1,km,2)=aa2(i,km)
             tauir(i,1,km)=bb(i,km)
  120      continue

c          do k=1,nadd
c            tausw(i,1,k)=0.0
c            tauir(i,1,k)=0.0
c          enddo
        enddo
      IF (IRADAVE .EQ. 1) THEN
         RNX1=1./FLOAT(NX1)
        DO K=1,NW
          Y1(K)=0.
          Y2(K)=0.
          Y3(K)=0.
          Y4(K)=0.
          Y5(K)=0.
          Y6(K)=0.
          Y7(K)=0.
        ENDDO
        DO 300 K=2,KLES
          KM=KMAX-K+NADD
          DO I=1,NX1
             Y1(KM)=Y1(KM)+TAQ(I,1,KM)
             Y2(KM)=Y2(KM)+WAQ(I,1,KM)
             Y3(KM)=Y3(KM)+TAUSW(I,1,KM,1)
             Y4(KM)=Y4(KM)+TAUSW(I,1,KM,2)
             Y5(KM)=Y5(KM)+REFF(I,1,KM,1)
             Y6(KM)=Y6(KM)+REFF(I,1,KM,2)
             Y7(KM)=Y7(KM)+TAUIR(I,1,KM)
          ENDDO
  300   CONTINUE
         DO K=2,KLES
          KM=KMAX-K+NADD
          TAQ(1,1,KM)=Y1(KM)*RNX1
          WAQ(1,1,KM)=Y2(KM)*RNX1
          TAUSW(1,1,KM,1)=Y3(KM)*RNX1
          TAUSW(1,1,KM,2)=Y4(KM)*RNX1
          REFF(1,1,KM,1)=Y5(KM)*RNX1
          REFF(1,1,KM,2)=Y6(KM)*RNX1
          TAUIR(1,1,KM)=Y7(KM)*RNX1
        ENDDO
      ENDIF
      return
      end

c********* u1/frmdc/radmod/cloud.f *****************                          
c-----lr,lt and lr2 are for an asymmetry factor 0f .843
      block data                                                                
      implicit none
      real    auc(26,21),c1(26,21),c2(26,21)
      common/co2/ auc,c1,c2
      real    auo(26,21),o1(26,21),o2(26,21)
      common/o3/ auo,o1,o2
      real    cah(22,19)
      common/sco2/ cah
      integer lr(10,10,31),lt(10,10,31),lr2(10,10)
      common/rtdata/lr,lt,lr2
      integer i,ip,iw,j,k
      data ((auc(ip,iw),iw=1,21),ip=1,3)/                                       
     &-0.3241e+01,-0.3099e+01,-0.2956e+01,-0.2810e+01,-0.2667e+01,              
     &-0.2532e+01,-0.2407e+01,-0.2290e+01,-0.2177e+01,-0.2067e+01,              
     &-0.1959e+01,-0.1852e+01,-0.1747e+01,-0.1643e+01,-0.1540e+01,              
     &-0.1439e+01,-0.1341e+01,-0.1247e+01,-0.1154e+01,-0.1063e+01,              
     &-0.9730e+00,-0.3241e+01,-0.3099e+01,-0.2955e+01,-0.2809e+01,              
     &-0.2666e+01,-0.2531e+01,-0.2405e+01,-0.2287e+01,-0.2173e+01,              
     &-0.2061e+01,-0.1951e+01,-0.1842e+01,-0.1735e+01,-0.1628e+01,              
     &-0.1523e+01,-0.1420e+01,-0.1320e+01,-0.1222e+01,-0.1127e+01,              
     &-0.1033e+01,-0.9400e+00,-0.3240e+01,-0.3098e+01,-0.2954e+01,              
     &-0.2808e+01,-0.2664e+01,-0.2528e+01,-0.2401e+01,-0.2282e+01,              
     &-0.2166e+01,-0.2052e+01,-0.1940e+01,-0.1828e+01,-0.1718e+01,              
     &-0.1609e+01,-0.1501e+01,-0.1396e+01,-0.1293e+01,-0.1193e+01,              
     &-0.1095e+01,-0.9980e+00,-0.9020e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=4,6)/                                       
     &-0.3239e+01,-0.3096e+01,-0.2952e+01,-0.2805e+01,-0.2661e+01,              
     &-0.2524e+01,-0.2396e+01,-0.2274e+01,-0.2156e+01,-0.2040e+01,              
     &-0.1924e+01,-0.1810e+01,-0.1698e+01,-0.1586e+01,-0.1475e+01,              
     &-0.1367e+01,-0.1262e+01,-0.1159e+01,-0.1057e+01,-0.9580e+00,              
     &-0.8600e+00,-0.3238e+01,-0.3094e+01,-0.2949e+01,-0.2802e+01,              
     &-0.2656e+01,-0.2518e+01,-0.2387e+01,-0.2263e+01,-0.2143e+01,              
     &-0.2023e+01,-0.1904e+01,-0.1787e+01,-0.1672e+01,-0.1557e+01,              
     &-0.1444e+01,-0.1334e+01,-0.1225e+01,-0.1119e+01,-0.1015e+01,              
     &-0.9130e+00,-0.8150e+00,-0.3235e+01,-0.3091e+01,-0.2945e+01,              
     &-0.2796e+01,-0.2649e+01,-0.2508e+01,-0.2375e+01,-0.2248e+01,              
     &-0.2124e+01,-0.2001e+01,-0.1879e+01,-0.1759e+01,-0.1640e+01,              
     &-0.1523e+01,-0.1407e+01,-0.1294e+01,-0.1183e+01,-0.1074e+01,              
     &-0.9680e+00,-0.8650e+00,-0.7670e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=7,9)/                                       
     &-0.3231e+01,-0.3085e+01,-0.2938e+01,-0.2788e+01,-0.2638e+01,              
     &-0.2495e+01,-0.2358e+01,-0.2227e+01,-0.2099e+01,-0.1972e+01,              
     &-0.1847e+01,-0.1724e+01,-0.1603e+01,-0.1483e+01,-0.1365e+01,              
     &-0.1249e+01,-0.1135e+01,-0.1024e+01,-0.9170e+00,-0.8150e+00,              
     &-0.7180e+00,-0.3225e+01,-0.3078e+01,-0.2928e+01,-0.2775e+01,              
     &-0.2623e+01,-0.2476e+01,-0.2335e+01,-0.2200e+01,-0.2068e+01,              
     &-0.1938e+01,-0.1809e+01,-0.1683e+01,-0.1559e+01,-0.1438e+01,              
     &-0.1317e+01,-0.1199e+01,-0.1084e+01,-0.9720e+00,-0.8660e+00,              
     &-0.7650e+00,-0.6730e+00,-0.3216e+01,-0.3066e+01,-0.2913e+01,              
     &-0.2757e+01,-0.2601e+01,-0.2450e+01,-0.2305e+01,-0.2165e+01,              
     &-0.2030e+01,-0.1896e+01,-0.1765e+01,-0.1636e+01,-0.1510e+01,              
     &-0.1386e+01,-0.1263e+01,-0.1144e+01,-0.1027e+01,-0.9160e+00,              
     &-0.8110e+00,-0.7130e+00,-0.6250e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=10,12)/                                     
     &-0.3202e+01,-0.3048e+01,-0.2892e+01,-0.2731e+01,-0.2571e+01,              
     &-0.2415e+01,-0.2266e+01,-0.2123e+01,-0.1984e+01,-0.1847e+01,              
     &-0.1713e+01,-0.1582e+01,-0.1454e+01,-0.1328e+01,-0.1205e+01,              
     &-0.1084e+01,-0.9680e+00,-0.8580e+00,-0.7550e+00,-0.6620e+00,              
     &-0.5790e+00,-0.3183e+01,-0.3023e+01,-0.2862e+01,-0.2697e+01,              
     &-0.2533e+01,-0.2372e+01,-0.2219e+01,-0.2073e+01,-0.1930e+01,              
     &-0.1791e+01,-0.1656e+01,-0.1523e+01,-0.1394e+01,-0.1266e+01,              
     &-0.1142e+01,-0.1021e+01,-0.9060e+00,-0.7990e+00,-0.7010e+00,              
     &-0.6130e+00,-0.5350e+00,-0.3156e+01,-0.2990e+01,-0.2823e+01,              
     &-0.2655e+01,-0.2486e+01,-0.2321e+01,-0.2165e+01,-0.2015e+01,              
     &-0.1871e+01,-0.1730e+01,-0.1593e+01,-0.1459e+01,-0.1329e+01,              
     &-0.1200e+01,-0.1076e+01,-0.9560e+00,-0.8440e+00,-0.7410e+00,              
     &-0.6480e+00,-0.5660e+00,-0.4940e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=13,15)/                                     
     &-0.3120e+01,-0.2948e+01,-0.2776e+01,-0.2603e+01,-0.2431e+01,              
     &-0.2263e+01,-0.2104e+01,-0.1952e+01,-0.1806e+01,-0.1664e+01,              
     &-0.1526e+01,-0.1392e+01,-0.1260e+01,-0.1132e+01,-0.1008e+01,              
     &-0.8920e+00,-0.7830e+00,-0.6850e+00,-0.5980e+00,-0.5220e+00,              
     &-0.4560e+00,-0.3078e+01,-0.2897e+01,-0.2720e+01,-0.2544e+01,              
     &-0.2369e+01,-0.2200e+01,-0.2038e+01,-0.1884e+01,-0.1737e+01,              
     &-0.1595e+01,-0.1456e+01,-0.1321e+01,-0.1189e+01,-0.1062e+01,              
     &-0.9410e+00,-0.8280e+00,-0.7250e+00,-0.6320e+00,-0.5520e+00,              
     &-0.4820e+00,-0.4210e+00,-0.3029e+01,-0.2840e+01,-0.2658e+01,              
     &-0.2480e+01,-0.2303e+01,-0.2132e+01,-0.1969e+01,-0.1814e+01,              
     &-0.1666e+01,-0.1524e+01,-0.1385e+01,-0.1249e+01,-0.1118e+01,              
     &-0.9930e+00,-0.8750e+00,-0.7660e+00,-0.6690e+00,-0.5830e+00,              
     &-0.5090e+00,-0.4450e+00,-0.3890e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=16,18)/                                     
     &-0.2978e+01,-0.2779e+01,-0.2593e+01,-0.2412e+01,-0.2235e+01,              
     &-0.2063e+01,-0.1899e+01,-0.1743e+01,-0.1595e+01,-0.1451e+01,              
     &-0.1311e+01,-0.1176e+01,-0.1046e+01,-0.9230e+00,-0.8090e+00,              
     &-0.7060e+00,-0.6140e+00,-0.5330e+00,-0.4630e+00,-0.4010e+00,              
     &-0.3480e+00,-0.2927e+01,-0.2717e+01,-0.2525e+01,-0.2343e+01,              
     &-0.2166e+01,-0.1995e+01,-0.1831e+01,-0.1674e+01,-0.1524e+01,              
     &-0.1379e+01,-0.1239e+01,-0.1104e+01,-0.9760e+00,-0.8560e+00,              
     &-0.7470e+00,-0.6490e+00,-0.5620e+00,-0.4850e+00,-0.4180e+00,              
     &-0.3580e+00,-0.3070e+00,-0.2883e+01,-0.2662e+01,-0.2461e+01,              
     &-0.2275e+01,-0.2098e+01,-0.1927e+01,-0.1762e+01,-0.1603e+01,              
     &-0.1450e+01,-0.1303e+01,-0.1161e+01,-0.1026e+01,-0.8980e+00,              
     &-0.7800e+00,-0.6720e+00,-0.5760e+00,-0.4930e+00,-0.4220e+00,              
     &-0.3620e+00,-0.3110e+00,-0.2660e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=19,21)/                                     
     &-0.2846e+01,-0.2612e+01,-0.2401e+01,-0.2210e+01,-0.2031e+01,              
     &-0.1861e+01,-0.1696e+01,-0.1535e+01,-0.1379e+01,-0.1228e+01,              
     &-0.1085e+01,-0.9500e+00,-0.8240e+00,-0.7090e+00,-0.6050e+00,              
     &-0.5160e+00,-0.4400e+00,-0.3770e+00,-0.3230e+00,-0.2760e+00,              
     &-0.2330e+00,-0.2817e+01,-0.2572e+01,-0.2350e+01,-0.2150e+01,              
     &-0.1968e+01,-0.1797e+01,-0.1632e+01,-0.1469e+01,-0.1310e+01,              
     &-0.1157e+01,-0.1013e+01,-0.8800e+00,-0.7580e+00,-0.6480e+00,              
     &-0.5520e+00,-0.4700e+00,-0.4020e+00,-0.3440e+00,-0.2930e+00,              
     &-0.2470e+00,-0.2060e+00,-0.2796e+01,-0.2541e+01,-0.2308e+01,              
     &-0.2098e+01,-0.1910e+01,-0.1736e+01,-0.1569e+01,-0.1405e+01,              
     &-0.1243e+01,-0.1089e+01,-0.9440e+00,-0.8120e+00,-0.6940e+00,              
     &-0.5910e+00,-0.5030e+00,-0.4300e+00,-0.3680e+00,-0.3140e+00,              
     &-0.2660e+00,-0.2220e+00,-0.1830e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=22,24)/                                     
     &-0.2781e+01,-0.2519e+01,-0.2275e+01,-0.2055e+01,-0.1857e+01,              
     &-0.1677e+01,-0.1507e+01,-0.1341e+01,-0.1179e+01,-0.1024e+01,              
     &-0.8800e+00,-0.7500e+00,-0.6370e+00,-0.5420e+00,-0.4620e+00,              
     &-0.3970e+00,-0.3400e+00,-0.2880e+00,-0.2410e+00,-0.2000e+00,              
     &-0.1650e+00,-0.2772e+01,-0.2504e+01,-0.2252e+01,-0.2022e+01,              
     &-0.1813e+01,-0.1624e+01,-0.1448e+01,-0.1280e+01,-0.1118e+01,              
     &-0.9630e+00,-0.8210e+00,-0.6950e+00,-0.5880e+00,-0.5000e+00,              
     &-0.4280e+00,-0.3680e+00,-0.3140e+00,-0.2640e+00,-0.2200e+00,              
     &-0.1810e+00,-0.1500e+00,-0.2766e+01,-0.2494e+01,-0.2236e+01,              
     &-0.1997e+01,-0.1778e+01,-0.1577e+01,-0.1393e+01,-0.1221e+01,              
     &-0.1058e+01,-0.9060e+00,-0.7670e+00,-0.6470e+00,-0.5460e+00,              
     &-0.4650e+00,-0.3990e+00,-0.3420e+00,-0.2900e+00,-0.2430e+00,              
     &-0.2010e+00,-0.1660e+00,-0.1380e+00/                                      
      data ((auc(ip,iw),iw=1,21),ip=25,26)/                                     
     &-0.2762e+01,-0.2486e+01,-0.2224e+01,-0.1977e+01,-0.1749e+01,              
     &-0.1538e+01,-0.1343e+01,-0.1164e+01,-0.1001e+01,-0.8520e+00,              
     &-0.7190e+00,-0.6050e+00,-0.5110e+00,-0.4360e+00,-0.3730e+00,              
     &-0.3180e+00,-0.2690e+00,-0.2250e+00,-0.1860e+00,-0.1550e+00,              
     &-0.1280e+00,-0.2759e+01,-0.2480e+01,-0.2213e+01,-0.1961e+01,              
     &-0.1725e+01,-0.1505e+01,-0.1299e+01,-0.1113e+01,-0.9470e+00,              
     &-0.8010e+00,-0.6750e+00,-0.5700e+00,-0.4850e+00,-0.4140e+00,              
     &-0.3520e+00,-0.2980e+00,-0.2500e+00,-0.2090e+00,-0.1750e+00,              
     &-0.1470e+00,-0.1210e+00/                                                  
c                                                                               
      data ((c1(ip,iw),iw=1,21),ip=1,3)/                                        
     & 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1200e-04, 0.1200e-04,              
     & 0.1300e-04, 0.1600e-04, 0.1900e-04, 0.2200e-04, 0.2300e-04,              
     & 0.2400e-04, 0.2500e-04, 0.2600e-04, 0.2700e-04, 0.2700e-04,              
     & 0.2500e-04, 0.2400e-04, 0.2200e-04, 0.2100e-04, 0.1900e-04,              
     & 0.1700e-04, 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1200e-04,              
     & 0.1200e-04, 0.1300e-04, 0.1600e-04, 0.1900e-04, 0.2200e-04,              
     & 0.2300e-04, 0.2300e-04, 0.2400e-04, 0.2500e-04, 0.2600e-04,              
     & 0.2600e-04, 0.2400e-04, 0.2200e-04, 0.2100e-04, 0.1900e-04,              
     & 0.1800e-04, 0.1600e-04, 0.6000e-05, 0.9000e-05, 0.1100e-04,              
     & 0.1200e-04, 0.1200e-04, 0.1300e-04, 0.1600e-04, 0.1900e-04,              
     & 0.2100e-04, 0.2200e-04, 0.2300e-04, 0.2300e-04, 0.2400e-04,              
     & 0.2500e-04, 0.2400e-04, 0.2300e-04, 0.2100e-04, 0.2000e-04,              
     & 0.1800e-04, 0.1600e-04, 0.1400e-04/                                      
      data ((c1(ip,iw),iw=1,21),ip=4,6)/                                        
     & 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1200e-04, 0.1200e-04,              
     & 0.1300e-04, 0.1600e-04, 0.1900e-04, 0.2100e-04, 0.2200e-04,              
     & 0.2200e-04, 0.2200e-04, 0.2300e-04, 0.2300e-04, 0.2300e-04,              
     & 0.2100e-04, 0.2000e-04, 0.1800e-04, 0.1700e-04, 0.1500e-04,              
     & 0.1300e-04, 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1200e-04,              
     & 0.1200e-04, 0.1300e-04, 0.1600e-04, 0.1800e-04, 0.2000e-04,              
     & 0.2100e-04, 0.2100e-04, 0.2100e-04, 0.2100e-04, 0.2200e-04,              
     & 0.2100e-04, 0.2000e-04, 0.1800e-04, 0.1700e-04, 0.1500e-04,              
     & 0.1400e-04, 0.1200e-04, 0.6000e-05, 0.9000e-05, 0.1100e-04,              
     & 0.1200e-04, 0.1200e-04, 0.1300e-04, 0.1500e-04, 0.1800e-04,              
     & 0.1900e-04, 0.2000e-04, 0.1900e-04, 0.1900e-04, 0.2000e-04,              
     & 0.2000e-04, 0.2000e-04, 0.1800e-04, 0.1700e-04, 0.1500e-04,              
     & 0.1400e-04, 0.1300e-04, 0.1100e-04/                                      
      data ((c1(ip,iw),iw=1,21),ip=7,9)/                                        
     & 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1200e-04, 0.1200e-04,              
     & 0.1300e-04, 0.1500e-04, 0.1700e-04, 0.1800e-04, 0.1800e-04,              
     & 0.1800e-04, 0.1800e-04, 0.1800e-04, 0.1800e-04, 0.1800e-04,              
     & 0.1700e-04, 0.1500e-04, 0.1400e-04, 0.1300e-04, 0.1100e-04,              
     & 0.1000e-04, 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1200e-04,              
     & 0.1200e-04, 0.1200e-04, 0.1400e-04, 0.1600e-04, 0.1700e-04,              
     & 0.1700e-04, 0.1600e-04, 0.1600e-04, 0.1700e-04, 0.1700e-04,              
     & 0.1700e-04, 0.1600e-04, 0.1500e-04, 0.1400e-04, 0.1200e-04,
     & 0.1100e-04, 0.1000e-04, 0.6000e-05, 0.9000e-05, 0.1100e-04,              
     & 0.1100e-04, 0.1100e-04, 0.1200e-04, 0.1300e-04, 0.1400e-04,              
     & 0.1500e-04, 0.1500e-04, 0.1500e-04, 0.1500e-04, 0.1500e-04,              
     & 0.1500e-04, 0.1500e-04, 0.1400e-04, 0.1300e-04, 0.1200e-04,              
     & 0.1100e-04, 0.1000e-04, 0.9000e-05/                                      
      data ((c1(ip,iw),iw=1,21),ip=10,12)/                                      
     & 0.6000e-05, 0.9000e-05, 0.1100e-04, 0.1100e-04, 0.1100e-04,              
     & 0.1100e-04, 0.1200e-04, 0.1300e-04, 0.1400e-04, 0.1400e-04,              
     & 0.1300e-04, 0.1300e-04, 0.1400e-04, 0.1400e-04, 0.1400e-04,              
     & 0.1300e-04, 0.1200e-04, 0.1200e-04, 0.1100e-04, 0.1000e-04,              
     & 0.8000e-05, 0.6000e-05, 0.8000e-05, 0.1000e-04, 0.1000e-04,              
     & 0.1000e-04, 0.1000e-04, 0.1000e-04, 0.1200e-04, 0.1200e-04,              
     & 0.1200e-04, 0.1200e-04, 0.1200e-04, 0.1300e-04, 0.1300e-04,              
     & 0.1300e-04, 0.1200e-04, 0.1200e-04, 0.1100e-04, 0.1000e-04,              
     & 0.9000e-05, 0.7000e-05, 0.6000e-05, 0.8000e-05, 0.9000e-05,              
     & 0.1000e-04, 0.9000e-05, 0.9000e-05, 0.9000e-05, 0.1000e-04,              
     & 0.1100e-04, 0.1100e-04, 0.1100e-04, 0.1100e-04, 0.1200e-04,              
     & 0.1200e-04, 0.1200e-04, 0.1100e-04, 0.1100e-04, 0.1000e-04,              
     & 0.9000e-05, 0.8000e-05, 0.6000e-05/                                      
      data ((c1(ip,iw),iw=1,21),ip=13,15)/                                      
     & 0.6000e-05, 0.8000e-05, 0.9000e-05, 0.9000e-05, 0.8000e-05,              
     & 0.7000e-05, 0.8000e-05, 0.9000e-05, 0.1000e-04, 0.1000e-04,              
     & 0.1000e-04, 0.1000e-04, 0.1100e-04, 0.1100e-04, 0.1100e-04,              
     & 0.1100e-04, 0.1000e-04, 0.9000e-05, 0.8000e-05, 0.6000e-05,              
     & 0.4000e-05, 0.5000e-05, 0.7000e-05, 0.8000e-05, 0.8000e-05,              
     & 0.7000e-05, 0.7000e-05, 0.7000e-05, 0.8000e-05, 0.9000e-05,              
     & 0.9000e-05, 0.1000e-04, 0.1000e-04, 0.1000e-04, 0.1000e-04,              
     & 0.1000e-04, 0.1000e-04, 0.1000e-04, 0.8000e-05, 0.7000e-05,              
     & 0.5000e-05, 0.3000e-05, 0.5000e-05, 0.6000e-05, 0.7000e-05,              
     & 0.7000e-05, 0.7000e-05, 0.6000e-05, 0.6000e-05, 0.7000e-05,              
     & 0.8000e-05, 0.9000e-05, 0.9000e-05, 0.1000e-04, 0.1000e-04,              
     & 0.1000e-04, 0.1000e-04, 0.1000e-04, 0.9000e-05, 0.7000e-05,              
     & 0.6000e-05, 0.4000e-05, 0.2000e-05/                                      
      data ((c1(ip,iw),iw=1,21),ip=16,18)/                                      
     & 0.4000e-05, 0.5000e-05, 0.6000e-05, 0.6000e-05, 0.6000e-05,              
     & 0.6000e-05, 0.5000e-05, 0.6000e-05, 0.7000e-05, 0.8000e-05,              
     & 0.8000e-05, 0.9000e-05, 0.9000e-05, 0.9000e-05, 0.9000e-05,              
     & 0.9000e-05, 0.8000e-05, 0.6000e-05, 0.4000e-05, 0.2000e-05,              
     & 0.1000e-05, 0.3000e-05, 0.4000e-05, 0.5000e-05, 0.6000e-05,              
     & 0.6000e-05, 0.6000e-05, 0.5000e-05, 0.5000e-05, 0.6000e-05,              
     & 0.7000e-05, 0.8000e-05, 0.8000e-05, 0.8000e-05, 0.9000e-05,              
     & 0.8000e-05, 0.8000e-05, 0.6000e-05, 0.4000e-05, 0.2000e-05,              
     & 0.1000e-05, 0.0000e+00, 0.2000e-05, 0.3000e-05, 0.4000e-05,              
     & 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.5000e-05,              
     & 0.5000e-05, 0.6000e-05, 0.6000e-05, 0.6000e-05, 0.6000e-05,              
     & 0.6000e-05, 0.6000e-05, 0.5000e-05, 0.4000e-05, 0.2000e-05,              
     & 0.1000e-05, 0.0000e+00, 0.0000e+00/                                      
      data ((c1(ip,iw),iw=1,21),ip=19,21)/                                      
     & 0.2000e-05, 0.2000e-05, 0.3000e-05, 0.4000e-05, 0.5000e-05,              
     & 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.4000e-05, 0.4000e-05,              
     & 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.4000e-05,              
     & 0.3000e-05, 0.1000e-05, 0.0000e+00, 0.0000e+00,-0.1000e-05,              
     &-0.2000e-05, 0.1000e-05, 0.2000e-05, 0.2000e-05, 0.3000e-05,              
     & 0.4000e-05, 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.4000e-05,              
     & 0.4000e-05, 0.4000e-05, 0.4000e-05, 0.5000e-05, 0.4000e-05,              
     & 0.3000e-05, 0.2000e-05, 0.0000e+00, 0.0000e+00,-0.1000e-05,              
     &-0.2000e-05,-0.3000e-05, 0.1000e-05, 0.1000e-05, 0.1000e-05,              
     & 0.2000e-05, 0.3000e-05, 0.4000e-05, 0.5000e-05, 0.5000e-05,              
     & 0.5000e-05, 0.4000e-05, 0.3000e-05, 0.3000e-05, 0.3000e-05,              
     & 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,-0.1000e-05,              
     &-0.2000e-05,-0.3000e-05,-0.3000e-05/                                      
      data ((c1(ip,iw),iw=1,21),ip=22,24)/                                      
     & 0.0000e+00, 0.0000e+00, 0.1000e-05, 0.1000e-05, 0.3000e-05,              
     & 0.4000e-05, 0.5000e-05, 0.5000e-05, 0.5000e-05, 0.4000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.2000e-05, 0.1000e-05,              
     & 0.0000e+00, 0.0000e+00,-0.1000e-05,-0.2000e-05,-0.3000e-05,              
     &-0.3000e-05, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.1000e-05,              
     & 0.2000e-05, 0.3000e-05, 0.4000e-05, 0.5000e-05, 0.5000e-05,              
     & 0.4000e-05, 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.1000e-05,              
     & 0.0000e+00, 0.0000e+00,-0.1000e-05,-0.2000e-05,-0.3000e-05,              
     &-0.3000e-05,-0.2000e-05, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.1000e-05, 0.2000e-05, 0.2000e-05, 0.3000e-05, 0.4000e-05,              
     & 0.5000e-05, 0.4000e-05, 0.3000e-05, 0.2000e-05, 0.1000e-05,              
     & 0.0000e+00,-0.1000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.3000e-05,-0.2000e-05,-0.1000e-05/                                      
      data ((c1(ip,iw),iw=1,21),ip=25,26)/                                      
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.1000e-05, 0.1000e-05,              
     & 0.2000e-05, 0.3000e-05, 0.4000e-05, 0.4000e-05, 0.4000e-05,              
     & 0.3000e-05, 0.2000e-05, 0.0000e+00,-0.1000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.1000e-05,              
     & 0.1000e-05, 0.2000e-05, 0.2000e-05, 0.3000e-05, 0.3000e-05,              
     & 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,-0.1000e-05,              
     &-0.3000e-05,-0.4000e-05,-0.3000e-05,-0.3000e-05,-0.1000e-05,              
     & 0.0000e+00, 0.0000e+00/                                                  
c                                                                               
      data ((c2(ip,iw),iw=1,21),ip=1,3)/                                        
     & 0.4300e-02, 0.5300e-02, 0.6000e-02, 0.6500e-02, 0.6700e-02,              
     & 0.6700e-02, 0.6800e-02, 0.7000e-02, 0.7400e-02, 0.7700e-02,              
     & 0.7900e-02, 0.8100e-02, 0.8200e-02, 0.8200e-02, 0.8200e-02,              
     & 0.8100e-02, 0.7900e-02, 0.7700e-02, 0.7400e-02, 0.7200e-02,              
     & 0.6800e-02, 0.4300e-02, 0.5300e-02, 0.6000e-02, 0.6500e-02,              
     & 0.6600e-02, 0.6700e-02, 0.6800e-02, 0.7000e-02, 0.7200e-02,              
     & 0.7500e-02, 0.7700e-02, 0.7900e-02, 0.7900e-02, 0.8000e-02,              
     & 0.7900e-02, 0.7700e-02, 0.7600e-02, 0.7300e-02, 0.7100e-02,              
     & 0.6800e-02, 0.6400e-02, 0.4300e-02, 0.5300e-02, 0.6000e-02,              
     & 0.6400e-02, 0.6600e-02, 0.6600e-02, 0.6700e-02, 0.6900e-02,              
     & 0.7100e-02, 0.7300e-02, 0.7500e-02, 0.7600e-02, 0.7600e-02,              
     & 0.7600e-02, 0.7600e-02, 0.7400e-02, 0.7200e-02, 0.7000e-02,              
     & 0.6700e-02, 0.6400e-02, 0.6000e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=4,6)/                                        
     & 0.4300e-02, 0.5200e-02, 0.6000e-02, 0.6400e-02, 0.6500e-02,              
     & 0.6500e-02, 0.6600e-02, 0.6700e-02, 0.6900e-02, 0.7100e-02,              
     & 0.7200e-02, 0.7300e-02, 0.7300e-02, 0.7300e-02, 0.7200e-02,              
     & 0.7000e-02, 0.6800e-02, 0.6600e-02, 0.6300e-02, 0.6000e-02,              
     & 0.5700e-02, 0.4200e-02, 0.5200e-02, 0.5900e-02, 0.6300e-02,              
     & 0.6400e-02, 0.6400e-02, 0.6400e-02, 0.6500e-02, 0.6600e-02,              
     & 0.6800e-02, 0.6900e-02, 0.6900e-02, 0.6900e-02, 0.6900e-02,              
     & 0.6800e-02, 0.6600e-02, 0.6400e-02, 0.6100e-02, 0.5900e-02,              
     & 0.5600e-02, 0.5300e-02, 0.4200e-02, 0.5100e-02, 0.5800e-02,              
     & 0.6200e-02, 0.6300e-02, 0.6200e-02, 0.6200e-02, 0.6200e-02,              
     & 0.6300e-02, 0.6400e-02, 0.6500e-02, 0.6500e-02, 0.6500e-02,              
     & 0.6500e-02, 0.6400e-02, 0.6200e-02, 0.6000e-02, 0.5700e-02,              
     & 0.5500e-02, 0.5200e-02, 0.5000e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=7,9)/                                        
     & 0.4100e-02, 0.5000e-02, 0.5700e-02, 0.6000e-02, 0.6000e-02,              
     & 0.5900e-02, 0.5900e-02, 0.5900e-02, 0.5900e-02, 0.6000e-02,              
     & 0.6100e-02, 0.6100e-02, 0.6100e-02, 0.6100e-02, 0.6000e-02,              
     & 0.5800e-02, 0.5600e-02, 0.5400e-02, 0.5100e-02, 0.4900e-02,              
     & 0.4700e-02, 0.4000e-02, 0.4800e-02, 0.5400e-02, 0.5700e-02,              
     & 0.5700e-02, 0.5600e-02, 0.5500e-02, 0.5500e-02, 0.5500e-02,              
     & 0.5600e-02, 0.5700e-02, 0.5700e-02, 0.5700e-02, 0.5600e-02,              
     & 0.5500e-02, 0.5400e-02, 0.5200e-02, 0.5000e-02, 0.4800e-02,              
     & 0.4600e-02, 0.4400e-02, 0.3800e-02, 0.4600e-02, 0.5100e-02,              
     & 0.5400e-02, 0.5400e-02, 0.5200e-02, 0.5100e-02, 0.5100e-02,              
     & 0.5100e-02, 0.5200e-02, 0.5300e-02, 0.5300e-02, 0.5300e-02,              
     & 0.5200e-02, 0.5200e-02, 0.5000e-02, 0.4800e-02, 0.4700e-02,              
     & 0.4500e-02, 0.4300e-02, 0.4200e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=10,12)/                                      
     & 0.3600e-02, 0.4300e-02, 0.4800e-02, 0.5000e-02, 0.4900e-02,              
     & 0.4800e-02, 0.4700e-02, 0.4600e-02, 0.4700e-02, 0.4800e-02,              
     & 0.4900e-02, 0.4900e-02, 0.4900e-02, 0.4900e-02, 0.4800e-02,              
     & 0.4700e-02, 0.4500e-02, 0.4400e-02, 0.4200e-02, 0.4100e-02,              
     & 0.4000e-02, 0.3300e-02, 0.3900e-02, 0.4300e-02, 0.4500e-02,              
     & 0.4500e-02, 0.4300e-02, 0.4200e-02, 0.4300e-02, 0.4300e-02,              
     & 0.4400e-02, 0.4500e-02, 0.4600e-02, 0.4600e-02, 0.4600e-02,              
     & 0.4500e-02, 0.4400e-02, 0.4300e-02, 0.4200e-02, 0.4000e-02,              
     & 0.4000e-02, 0.3900e-02, 0.2900e-02, 0.3400e-02, 0.3800e-02,              
     & 0.4000e-02, 0.4000e-02, 0.3900e-02, 0.3900e-02, 0.3900e-02,              
     & 0.4000e-02, 0.4100e-02, 0.4200e-02, 0.4300e-02, 0.4300e-02,              
     & 0.4300e-02, 0.4200e-02, 0.4200e-02, 0.4100e-02, 0.4000e-02,              
     & 0.3900e-02, 0.3800e-02, 0.3800e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=13,15)/                                      
     & 0.2400e-02, 0.2900e-02, 0.3200e-02, 0.3500e-02, 0.3600e-02,              
     & 0.3600e-02, 0.3600e-02, 0.3600e-02, 0.3700e-02, 0.3900e-02,              
     & 0.4000e-02, 0.4000e-02, 0.4000e-02, 0.4000e-02, 0.4000e-02,              
     & 0.4000e-02, 0.3900e-02, 0.3800e-02, 0.3800e-02, 0.3700e-02,              
     & 0.3700e-02, 0.1900e-02, 0.2400e-02, 0.2700e-02, 0.3000e-02,              
     & 0.3200e-02, 0.3300e-02, 0.3300e-02, 0.3400e-02, 0.3500e-02,              
     & 0.3600e-02, 0.3700e-02, 0.3800e-02, 0.3800e-02, 0.3900e-02,              
     & 0.3800e-02, 0.3800e-02, 0.3800e-02, 0.3700e-02, 0.3700e-02,              
     & 0.3700e-02, 0.3600e-02, 0.1500e-02, 0.1900e-02, 0.2300e-02,              
     & 0.2600e-02, 0.2900e-02, 0.3000e-02, 0.3100e-02, 0.3200e-02,              
     & 0.3300e-02, 0.3500e-02, 0.3600e-02, 0.3600e-02, 0.3700e-02,              
     & 0.3700e-02, 0.3700e-02, 0.3700e-02, 0.3700e-02, 0.3600e-02,              
     & 0.3600e-02, 0.3600e-02, 0.3600e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=16,18)/                                      
     & 0.1100e-02, 0.1500e-02, 0.1900e-02, 0.2300e-02, 0.2600e-02,              
     & 0.2900e-02, 0.3000e-02, 0.3100e-02, 0.3200e-02, 0.3300e-02,              
     & 0.3400e-02, 0.3500e-02, 0.3600e-02, 0.3600e-02, 0.3600e-02,              
     & 0.3700e-02, 0.3700e-02, 0.3700e-02, 0.3700e-02, 0.3600e-02,              
     & 0.3500e-02, 0.8000e-03, 0.1100e-02, 0.1500e-02, 0.2000e-02,              
     & 0.2400e-02, 0.2700e-02, 0.2900e-02, 0.3000e-02, 0.3100e-02,              
     & 0.3200e-02, 0.3300e-02, 0.3400e-02, 0.3500e-02, 0.3600e-02,              
     & 0.3600e-02, 0.3600e-02, 0.3700e-02, 0.3700e-02, 0.3700e-02,              
     & 0.3600e-02, 0.3400e-02, 0.6000e-03, 0.9000e-03, 0.1300e-02,              
     & 0.1700e-02, 0.2100e-02, 0.2500e-02, 0.2800e-02, 0.2900e-02,              
     & 0.3000e-02, 0.3100e-02, 0.3200e-02, 0.3200e-02, 0.3300e-02,              
     & 0.3300e-02, 0.3300e-02, 0.3300e-02, 0.3300e-02, 0.3300e-02,              
     & 0.3300e-02, 0.3200e-02, 0.3100e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=19,21)/                                      
     & 0.4000e-03, 0.7000e-03, 0.1000e-02, 0.1400e-02, 0.1900e-02,              
     & 0.2300e-02, 0.2600e-02, 0.2800e-02, 0.2900e-02, 0.3000e-02,              
     & 0.3000e-02, 0.3100e-02, 0.3100e-02, 0.3200e-02, 0.3200e-02,              
     & 0.3200e-02, 0.3300e-02, 0.3300e-02, 0.3300e-02, 0.3200e-02,              
     & 0.3000e-02, 0.3000e-03, 0.5000e-03, 0.8000e-03, 0.1200e-02,              
     & 0.1700e-02, 0.2100e-02, 0.2400e-02, 0.2700e-02, 0.2800e-02,              
     & 0.2900e-02, 0.2900e-02, 0.3000e-02, 0.3000e-02, 0.3100e-02,              
     & 0.3200e-02, 0.3200e-02, 0.3300e-02, 0.3300e-02, 0.3200e-02,              
     & 0.3100e-02, 0.3000e-02, 0.2000e-03, 0.4000e-03, 0.7000e-03,              
     & 0.1100e-02, 0.1400e-02, 0.1800e-02, 0.2200e-02, 0.2500e-02,              
     & 0.2700e-02, 0.2800e-02, 0.2800e-02, 0.2800e-02, 0.2900e-02,              
     & 0.3000e-02, 0.3100e-02, 0.3200e-02, 0.3200e-02, 0.3200e-02,              
     & 0.3200e-02, 0.3000e-02, 0.2900e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=22,24)/                                      
     & 0.2000e-03, 0.4000e-03, 0.6000e-03, 0.9000e-03, 0.1200e-02,              
     & 0.1600e-02, 0.1900e-02, 0.2200e-02, 0.2500e-02, 0.2600e-02,              
     & 0.2700e-02, 0.2700e-02, 0.2800e-02, 0.2900e-02, 0.3000e-02,              
     & 0.3200e-02, 0.3200e-02, 0.3200e-02, 0.3000e-02, 0.2900e-02,              
     & 0.2700e-02, 0.2000e-03, 0.3000e-03, 0.6000e-03, 0.8000e-03,              
     & 0.1100e-02, 0.1300e-02, 0.1600e-02, 0.2000e-02, 0.2300e-02,              
     & 0.2500e-02, 0.2600e-02, 0.2700e-02, 0.2800e-02, 0.2900e-02,              
     & 0.3000e-02, 0.3100e-02, 0.3100e-02, 0.3100e-02, 0.2900e-02,              
     & 0.2700e-02, 0.2600e-02, 0.2000e-03, 0.3000e-03, 0.5000e-03,              
     & 0.7000e-03, 0.9000e-03, 0.1100e-02, 0.1400e-02, 0.1700e-02,              
     & 0.2000e-02, 0.2300e-02, 0.2500e-02, 0.2700e-02, 0.2800e-02,              
     & 0.2900e-02, 0.3000e-02, 0.3000e-02, 0.3000e-02, 0.2900e-02,              
     & 0.2800e-02, 0.2600e-02, 0.2500e-02/                                      
      data ((c2(ip,iw),iw=1,21),ip=25,26)/                                      
     & 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03, 0.8000e-03,              
     & 0.9000e-03, 0.1100e-02, 0.1400e-02, 0.1700e-02, 0.2000e-02,              
     & 0.2300e-02, 0.2600e-02, 0.2800e-02, 0.3000e-02, 0.3000e-02,              
     & 0.3000e-02, 0.2900e-02, 0.2800e-02, 0.2600e-02, 0.2500e-02,              
     & 0.2500e-02, 0.1000e-03, 0.2000e-03, 0.3000e-03, 0.5000e-03,              
     & 0.6000e-03, 0.7000e-03, 0.9000e-03, 0.1100e-02, 0.1500e-02,              
     & 0.1800e-02, 0.2200e-02, 0.2600e-02, 0.2900e-02, 0.3100e-02,              
     & 0.3000e-02, 0.2900e-02, 0.2800e-02, 0.2600e-02, 0.2500e-02,              
     & 0.2500e-02, 0.2500e-02/                                                  
c                                                                               
      data ((auo(ip,iw),iw=1,21),ip=1,3)/                                       
     &-0.3186e+01,-0.2996e+01,-0.2811e+01,-0.2632e+01,-0.2460e+01,              
     &-0.2297e+01,-0.2145e+01,-0.2007e+01,-0.1881e+01,-0.1768e+01,              
     &-0.1665e+01,-0.1569e+01,-0.1478e+01,-0.1389e+01,-0.1302e+01,              
     &-0.1216e+01,-0.1132e+01,-0.1051e+01,-0.9740e+00,-0.9010e+00,              
     &-0.8340e+00,-0.3186e+01,-0.2996e+01,-0.2811e+01,-0.2632e+01,              
     &-0.2460e+01,-0.2297e+01,-0.2145e+01,-0.2007e+01,-0.1881e+01,              
     &-0.1768e+01,-0.1665e+01,-0.1569e+01,-0.1478e+01,-0.1389e+01,              
     &-0.1301e+01,-0.1215e+01,-0.1131e+01,-0.1050e+01,-0.9720e+00,              
     &-0.8990e+00,-0.8320e+00,-0.3186e+01,-0.2996e+01,-0.2811e+01,              
     &-0.2632e+01,-0.2460e+01,-0.2297e+01,-0.2145e+01,-0.2006e+01,              
     &-0.1881e+01,-0.1768e+01,-0.1665e+01,-0.1568e+01,-0.1477e+01,              
     &-0.1388e+01,-0.1300e+01,-0.1214e+01,-0.1130e+01,-0.1048e+01,              
     &-0.9700e+00,-0.8970e+00,-0.8290e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=4,6)/                                       
     &-0.3186e+01,-0.2996e+01,-0.2811e+01,-0.2632e+01,-0.2459e+01,              
     &-0.2297e+01,-0.2145e+01,-0.2006e+01,-0.1880e+01,-0.1767e+01,              
     &-0.1664e+01,-0.1568e+01,-0.1476e+01,-0.1387e+01,-0.1299e+01,              
     &-0.1212e+01,-0.1127e+01,-0.1045e+01,-0.9670e+00,-0.8930e+00,              
     &-0.8240e+00,-0.3186e+01,-0.2996e+01,-0.2811e+01,-0.2631e+01,              
     &-0.2459e+01,-0.2296e+01,-0.2145e+01,-0.2005e+01,-0.1880e+01,              
     &-0.1766e+01,-0.1663e+01,-0.1566e+01,-0.1474e+01,-0.1385e+01,              
     &-0.1296e+01,-0.1209e+01,-0.1124e+01,-0.1041e+01,-0.9620e+00,              
     &-0.8870e+00,-0.8160e+00,-0.3186e+01,-0.2996e+01,-0.2811e+01,              
     &-0.2631e+01,-0.2459e+01,-0.2296e+01,-0.2144e+01,-0.2005e+01,              
     &-0.1879e+01,-0.1765e+01,-0.1661e+01,-0.1564e+01,-0.1471e+01,              
     &-0.1381e+01,-0.1292e+01,-0.1205e+01,-0.1118e+01,-0.1035e+01,              
     &-0.9540e+00,-0.8780e+00,-0.8060e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=7,9)/                                       
     &-0.3186e+01,-0.2996e+01,-0.2811e+01,-0.2631e+01,-0.2458e+01,              
     &-0.2295e+01,-0.2143e+01,-0.2003e+01,-0.1877e+01,-0.1763e+01,              
     &-0.1658e+01,-0.1561e+01,-0.1467e+01,-0.1376e+01,-0.1287e+01,              
     &-0.1198e+01,-0.1110e+01,-0.1025e+01,-0.9430e+00,-0.8640e+00,              
     &-0.7900e+00,-0.3185e+01,-0.2996e+01,-0.2810e+01,-0.2630e+01,              
     &-0.2458e+01,-0.2294e+01,-0.2141e+01,-0.2001e+01,-0.1874e+01,              
     &-0.1759e+01,-0.1654e+01,-0.1555e+01,-0.1461e+01,-0.1369e+01,              
     &-0.1278e+01,-0.1187e+01,-0.1098e+01,-0.1011e+01,-0.9270e+00,              
     &-0.8460e+00,-0.7700e+00,-0.3185e+01,-0.2996e+01,-0.2810e+01,              
     &-0.2629e+01,-0.2456e+01,-0.2292e+01,-0.2139e+01,-0.1998e+01,              
     &-0.1870e+01,-0.1754e+01,-0.1647e+01,-0.1547e+01,-0.1451e+01,              
     &-0.1357e+01,-0.1264e+01,-0.1172e+01,-0.1080e+01,-0.9910e+00,              
     &-0.9040e+00,-0.8210e+00,-0.7430e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=10,12)/                                     
     &-0.3185e+01,-0.2995e+01,-0.2809e+01,-0.2628e+01,-0.2454e+01,              
     &-0.2289e+01,-0.2135e+01,-0.1993e+01,-0.1863e+01,-0.1746e+01,              
     &-0.1637e+01,-0.1535e+01,-0.1437e+01,-0.1341e+01,-0.1245e+01,              
     &-0.1150e+01,-0.1056e+01,-0.9640e+00,-0.8750e+00,-0.7900e+00,              
     &-0.7090e+00,-0.3185e+01,-0.2994e+01,-0.2807e+01,-0.2626e+01,              
     &-0.2451e+01,-0.2285e+01,-0.2129e+01,-0.1985e+01,-0.1854e+01,              
     &-0.1733e+01,-0.1622e+01,-0.1518e+01,-0.1417e+01,-0.1317e+01,              
     &-0.1219e+01,-0.1121e+01,-0.1024e+01,-0.9300e+00,-0.8380e+00,              
     &-0.7510e+00,-0.6690e+00,-0.3183e+01,-0.2993e+01,-0.2805e+01,              
     &-0.2623e+01,-0.2447e+01,-0.2279e+01,-0.2121e+01,-0.1974e+01,              
     &-0.1840e+01,-0.1716e+01,-0.1602e+01,-0.1493e+01,-0.1389e+01,              
     &-0.1286e+01,-0.1184e+01,-0.1083e+01,-0.9840e+00,-0.8880e+00,              
     &-0.7950e+00,-0.7060e+00,-0.6230e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=13,15)/                                     
     &-0.3183e+01,-0.2991e+01,-0.2802e+01,-0.2619e+01,-0.2441e+01,              
     &-0.2270e+01,-0.2110e+01,-0.1959e+01,-0.1821e+01,-0.1693e+01,              
     &-0.1574e+01,-0.1461e+01,-0.1352e+01,-0.1246e+01,-0.1141e+01,              
     &-0.1038e+01,-0.9370e+00,-0.8390e+00,-0.7450e+00,-0.6560e+00,              
     &-0.5720e+00,-0.3181e+01,-0.2988e+01,-0.2799e+01,-0.2613e+01,              
     &-0.2433e+01,-0.2260e+01,-0.2095e+01,-0.1940e+01,-0.1796e+01,              
     &-0.1662e+01,-0.1538e+01,-0.1420e+01,-0.1307e+01,-0.1198e+01,              
     &-0.1090e+01,-0.9850e+00,-0.8830e+00,-0.7840e+00,-0.6900e+00,              
     &-0.6020e+00,-0.5200e+00,-0.3179e+01,-0.2985e+01,-0.2794e+01,              
     &-0.2606e+01,-0.2424e+01,-0.2247e+01,-0.2077e+01,-0.1917e+01,              
     &-0.1767e+01,-0.1627e+01,-0.1496e+01,-0.1373e+01,-0.1256e+01,              
     &-0.1143e+01,-0.1033e+01,-0.9270e+00,-0.8240e+00,-0.7260e+00,              
     &-0.6330e+00,-0.5460e+00,-0.4670e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=16,18)/                                     
     &-0.3177e+01,-0.2982e+01,-0.2789e+01,-0.2599e+01,-0.2413e+01,              
     &-0.2232e+01,-0.2058e+01,-0.1892e+01,-0.1734e+01,-0.1587e+01,              
     &-0.1450e+01,-0.1322e+01,-0.1200e+01,-0.1084e+01,-0.9730e+00,              
     &-0.8670e+00,-0.7640e+00,-0.6670e+00,-0.5760e+00,-0.4920e+00,              
     &-0.4160e+00,-0.3167e+01,-0.2971e+01,-0.2777e+01,-0.2585e+01,              
     &-0.2396e+01,-0.2211e+01,-0.2033e+01,-0.1861e+01,-0.1698e+01,              
     &-0.1544e+01,-0.1400e+01,-0.1267e+01,-0.1141e+01,-0.1023e+01,              
     &-0.9110e+00,-0.8040e+00,-0.7020e+00,-0.6070e+00,-0.5180e+00,              
     &-0.4370e+00,-0.3640e+00,-0.3169e+01,-0.2972e+01,-0.2776e+01,              
     &-0.2582e+01,-0.2390e+01,-0.2203e+01,-0.2019e+01,-0.1841e+01,              
     &-0.1671e+01,-0.1510e+01,-0.1358e+01,-0.1217e+01,-0.1086e+01,              
     &-0.9640e+00,-0.8500e+00,-0.7430e+00,-0.6430e+00,-0.5510e+00,              
     &-0.4650e+00,-0.3890e+00,-0.3210e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=19,21)/                                     
     &-0.3170e+01,-0.2972e+01,-0.2775e+01,-0.2580e+01,-0.2386e+01,              
     &-0.2195e+01,-0.2008e+01,-0.1826e+01,-0.1650e+01,-0.1481e+01,              
     &-0.1322e+01,-0.1174e+01,-0.1036e+01,-0.9100e+00,-0.7940e+00,              
     &-0.6870e+00,-0.5890e+00,-0.4990e+00,-0.4180e+00,-0.3460e+00,              
     &-0.2840e+00,-0.3171e+01,-0.2972e+01,-0.2774e+01,-0.2578e+01,              
     &-0.2382e+01,-0.2189e+01,-0.1999e+01,-0.1813e+01,-0.1633e+01,              
     &-0.1459e+01,-0.1293e+01,-0.1138e+01,-0.9940e+00,-0.8620e+00,              
     &-0.7430e+00,-0.6360e+00,-0.5400e+00,-0.4540e+00,-0.3770e+00,              
     &-0.3110e+00,-0.2530e+00,-0.3170e+01,-0.2971e+01,-0.2773e+01,              
     &-0.2576e+01,-0.2379e+01,-0.2185e+01,-0.1993e+01,-0.1804e+01,              
     &-0.1620e+01,-0.1441e+01,-0.1271e+01,-0.1110e+01,-0.9600e+00,              
     &-0.8230e+00,-0.7010e+00,-0.5920e+00,-0.4970e+00,-0.4140e+00,              
     &-0.3430e+00,-0.2810e+00,-0.2290e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=22,24)/                                     
     &-0.3170e+01,-0.2971e+01,-0.2772e+01,-0.2574e+01,-0.2377e+01,              
     &-0.2182e+01,-0.1988e+01,-0.1798e+01,-0.1611e+01,-0.1430e+01,              
     &-0.1255e+01,-0.1089e+01,-0.9340e+00,-0.7930e+00,-0.6670e+00,              
     &-0.5560e+00,-0.4620e+00,-0.3820e+00,-0.3150e+00,-0.2580e+00,              
     &-0.2110e+00,-0.3170e+01,-0.2971e+01,-0.2772e+01,-0.2573e+01,              
     &-0.2376e+01,-0.2180e+01,-0.1985e+01,-0.1793e+01,-0.1605e+01,              
     &-0.1421e+01,-0.1244e+01,-0.1075e+01,-0.9160e+00,-0.7710e+00,              
     &-0.6410e+00,-0.5290e+00,-0.4340e+00,-0.3560e+00,-0.2930e+00,              
     &-0.2410e+00,-0.1970e+00,-0.3170e+01,-0.2971e+01,-0.2771e+01,              
     &-0.2573e+01,-0.2375e+01,-0.2178e+01,-0.1983e+01,-0.1791e+01,              
     &-0.1601e+01,-0.1416e+01,-0.1237e+01,-0.1065e+01,-0.9040e+00,              
     &-0.7560e+00,-0.6230e+00,-0.5090e+00,-0.4140e+00,-0.3380e+00,              
     &-0.2780e+00,-0.2290e+00,-0.1880e+00/                                      
      data ((auo(ip,iw),iw=1,21),ip=25,26)/                                     
     &-0.3171e+01,-0.2971e+01,-0.2772e+01,-0.2573e+01,-0.2375e+01,              
     &-0.2178e+01,-0.1983e+01,-0.1790e+01,-0.1600e+01,-0.1414e+01,              
     &-0.1233e+01,-0.1060e+01,-0.8970e+00,-0.7470e+00,-0.6120e+00,              
     &-0.4960e+00,-0.4010e+00,-0.3260e+00,-0.2680e+00,-0.2220e+00,              
     &-0.1830e+00,-0.3172e+01,-0.2973e+01,-0.2773e+01,-0.2575e+01,              
     &-0.2376e+01,-0.2179e+01,-0.1984e+01,-0.1790e+01,-0.1600e+01,              
     &-0.1413e+01,-0.1232e+01,-0.1058e+01,-0.8930e+00,-0.7420e+00,              
     &-0.6060e+00,-0.4890e+00,-0.3940e+00,-0.3190e+00,-0.2620e+00,              
     &-0.2170e+00,-0.1790e+00/                                                  
c                                                                               
      data ((o1(ip,iw),iw=1,21),ip=1,3)/                                        
     &-0.3000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,-0.5000e-05,              
     &-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.4000e-05,-0.3000e-05,              
     &-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05, 0.3000e-05,              
     & 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,-0.1000e-05,              
     &-0.3000e-05,-0.3000e-05,-0.3000e-05,-0.3000e-05,-0.4000e-05,              
     &-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.4000e-05,              
     &-0.3000e-05,-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,              
     &-0.1000e-05,-0.3000e-05,-0.3000e-05,-0.3000e-05,-0.4000e-05,              
     &-0.4000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,              
     &-0.4000e-05,-0.3000e-05,-0.1000e-05, 0.0000e+00, 0.1000e-05,              
     & 0.2000e-05, 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.1000e-05,              
     & 0.0000e+00,-0.1000e-05,-0.3000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=4,6)/                                        
     &-0.3000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,-0.5000e-05,              
     &-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.4000e-05,-0.3000e-05,              
     &-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05, 0.3000e-05,              
     & 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,-0.1000e-05,              
     &-0.3000e-05,-0.3000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,              
     &-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.4000e-05,              
     &-0.3000e-05,-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,              
     &-0.1000e-05,-0.3000e-05,-0.2000e-05,-0.3000e-05,-0.4000e-05,              
     &-0.4000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,              
     &-0.4000e-05,-0.3000e-05,-0.1000e-05, 0.0000e+00, 0.1000e-05,              
     & 0.2000e-05, 0.3000e-05, 0.3000e-05, 0.3000e-05, 0.1000e-05,              
     & 0.0000e+00,-0.1000e-05,-0.2000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=7,9)/                                        
     &-0.2000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,-0.5000e-05,              
     &-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.4000e-05,-0.3000e-05,              
     &-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05, 0.3000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.0000e+00,-0.1000e-05,              
     &-0.2000e-05,-0.3000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,              
     &-0.4000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,-0.4000e-05,              
     &-0.3000e-05,-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.0000e+00,              
     &-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.3000e-05,-0.3000e-05,              
     &-0.4000e-05,-0.4000e-05,-0.5000e-05,-0.5000e-05,-0.5000e-05,              
     &-0.4000e-05,-0.2000e-05,-0.1000e-05, 0.0000e+00, 0.1000e-05,              
     & 0.2000e-05, 0.3000e-05, 0.3000e-05, 0.3000e-05, 0.2000e-05,              
     & 0.0000e+00,-0.1000e-05,-0.2000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=10,12)/                                      
     &-0.3000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,-0.4000e-05,              
     &-0.4000e-05,-0.5000e-05,-0.4000e-05,-0.3000e-05,-0.2000e-05,              
     &-0.1000e-05, 0.0000e+00, 0.1000e-05, 0.2000e-05, 0.3000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.1000e-05, 0.0000e+00,-0.1000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.3000e-05,-0.3000e-05,-0.3000e-05,              
     &-0.4000e-05,-0.4000e-05,-0.4000e-05,-0.4000e-05,-0.3000e-05,              
     &-0.2000e-05, 0.0000e+00, 0.0000e+00, 0.1000e-05, 0.2000e-05,              
     & 0.3000e-05, 0.3000e-05, 0.2000e-05, 0.1000e-05, 0.0000e+00,              
     &-0.1000e-05,-0.2000e-05,-0.3000e-05,-0.3000e-05,-0.3000e-05,              
     &-0.3000e-05,-0.3000e-05,-0.4000e-05,-0.4000e-05,-0.3000e-05,              
     &-0.3000e-05,-0.2000e-05, 0.0000e+00, 0.0000e+00, 0.1000e-05,              
     & 0.2000e-05, 0.2000e-05, 0.2000e-05, 0.2000e-05, 0.1000e-05,              
     & 0.0000e+00,-0.1000e-05,-0.3000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=13,15)/                                      
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.3000e-05,-0.3000e-05,              
     &-0.3000e-05,-0.3000e-05,-0.3000e-05,-0.2000e-05,-0.1000e-05,              
     & 0.0000e+00, 0.0000e+00, 0.1000e-05, 0.2000e-05, 0.2000e-05,              
     & 0.2000e-05, 0.1000e-05, 0.0000e+00, 0.0000e+00,-0.1000e-05,              
     &-0.3000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.3000e-05,-0.3000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.1000e-05, 0.0000e+00, 0.0000e+00, 0.1000e-05, 0.1000e-05,              
     & 0.1000e-05, 0.1000e-05, 0.1000e-05, 0.0000e+00, 0.0000e+00,              
     &-0.2000e-05,-0.3000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.1000e-05, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.1000e-05, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     &-0.1000e-05,-0.2000e-05,-0.3000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=16,18)/                                      
     &-0.1000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.1000e-05,              
     &-0.1000e-05, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00,-0.1000e-05,-0.2000e-05,              
     &-0.3000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.1000e-05,-0.1000e-05, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,-0.1000e-05,              
     &-0.2000e-05,-0.3000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     &-0.1000e-05,-0.2000e-05,-0.2000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=19,21)/                                      
     &-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.1000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,              
     & 0.0000e+00, 0.0000e+00,-0.1000e-05,-0.1000e-05,-0.1000e-05,              
     &-0.2000e-05,-0.1000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.1000e-05,              
     &-0.1000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,              
     &-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,              
     &-0.2000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,              
     &-0.1000e-05,-0.2000e-05,-0.2000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=22,24)/                                      
     &-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.2000e-05,              
     &-0.1000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.1000e-05,-0.1000e-05,-0.1000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.1000e-05,-0.2000e-05,-0.1000e-05,-0.1000e-05,              
     &-0.2000e-05,-0.1000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.3000e-05/                                      
      data ((o1(ip,iw),iw=1,21),ip=25,26)/                                      
     &-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.2000e-05,              
     &-0.1000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.3000e-05,-0.2000e-05,-0.1000e-05,-0.2000e-05,-0.1000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.1000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,-0.2000e-05,              
     &-0.2000e-05,-0.3000e-05/                                                  
c                                                                               
      data ((o2(ip,iw),iw=1,21),ip=1,3)/                                        
     & 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03, 0.9000e-03,              
     & 0.1300e-02, 0.1700e-02, 0.2200e-02, 0.2800e-02, 0.3300e-02,              
     & 0.3700e-02, 0.4200e-02, 0.4500e-02, 0.4900e-02, 0.5100e-02,              
     & 0.5400e-02, 0.5500e-02, 0.5500e-02, 0.5500e-02, 0.5300e-02,              
     & 0.5200e-02, 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.9000e-03, 0.1300e-02, 0.1700e-02, 0.2200e-02, 0.2800e-02,              
     & 0.3300e-02, 0.3700e-02, 0.4200e-02, 0.4500e-02, 0.4900e-02,              
     & 0.5100e-02, 0.5300e-02, 0.5500e-02, 0.5500e-02, 0.5400e-02,              
     & 0.5300e-02, 0.5100e-02, 0.1000e-03, 0.2000e-03, 0.4000e-03,              
     & 0.6000e-03, 0.9000e-03, 0.1300e-02, 0.1700e-02, 0.2200e-02,              
     & 0.2700e-02, 0.3300e-02, 0.3700e-02, 0.4100e-02, 0.4500e-02,              
     & 0.4800e-02, 0.5100e-02, 0.5300e-02, 0.5400e-02, 0.5400e-02,              
     & 0.5400e-02, 0.5200e-02, 0.5100e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=4,6)/                                        
     & 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03, 0.9000e-03,              
     & 0.1300e-02, 0.1700e-02, 0.2200e-02, 0.2700e-02, 0.3200e-02,              
     & 0.3700e-02, 0.4100e-02, 0.4500e-02, 0.4800e-02, 0.5100e-02,              
     & 0.5300e-02, 0.5400e-02, 0.5400e-02, 0.5300e-02, 0.5200e-02,              
     & 0.5000e-02, 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.9000e-03, 0.1300e-02, 0.1700e-02, 0.2200e-02, 0.2700e-02,              
     & 0.3200e-02, 0.3700e-02, 0.4100e-02, 0.4500e-02, 0.4800e-02,              
     & 0.5000e-02, 0.5200e-02, 0.5300e-02, 0.5300e-02, 0.5200e-02,              
     & 0.5000e-02, 0.4800e-02, 0.1000e-03, 0.2000e-03, 0.4000e-03,              
     & 0.6000e-03, 0.9000e-03, 0.1300e-02, 0.1700e-02, 0.2200e-02,              
     & 0.2700e-02, 0.3200e-02, 0.3700e-02, 0.4100e-02, 0.4400e-02,              
     & 0.4700e-02, 0.5000e-02, 0.5100e-02, 0.5200e-02, 0.5200e-02,              
     & 0.5100e-02, 0.4900e-02, 0.4600e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=7,9)/                                        
     & 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03, 0.9000e-03,              
     & 0.1300e-02, 0.1700e-02, 0.2200e-02, 0.2700e-02, 0.3200e-02,              
     & 0.3600e-02, 0.4000e-02, 0.4300e-02, 0.4600e-02, 0.4800e-02,              
     & 0.5000e-02, 0.5000e-02, 0.5000e-02, 0.4800e-02, 0.4600e-02,              
     & 0.4400e-02, 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.9000e-03, 0.1200e-02, 0.1600e-02, 0.2100e-02, 0.2600e-02,              
     & 0.3100e-02, 0.3500e-02, 0.3900e-02, 0.4200e-02, 0.4500e-02,              
     & 0.4700e-02, 0.4800e-02, 0.4800e-02, 0.4700e-02, 0.4600e-02,              
     & 0.4400e-02, 0.4100e-02, 0.1000e-03, 0.2000e-03, 0.4000e-03,              
     & 0.6000e-03, 0.9000e-03, 0.1200e-02, 0.1600e-02, 0.2100e-02,              
     & 0.2500e-02, 0.3000e-02, 0.3400e-02, 0.3800e-02, 0.4100e-02,              
     & 0.4300e-02, 0.4500e-02, 0.4500e-02, 0.4500e-02, 0.4400e-02,              
     & 0.4200e-02, 0.4000e-02, 0.3800e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=10,12)/                                      
     & 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03, 0.8000e-03,              
     & 0.1100e-02, 0.1500e-02, 0.2000e-02, 0.2400e-02, 0.2900e-02,              
     & 0.3200e-02, 0.3600e-02, 0.3800e-02, 0.4000e-02, 0.4100e-02,              
     & 0.4200e-02, 0.4200e-02, 0.4000e-02, 0.3900e-02, 0.3700e-02,              
     & 0.3400e-02, 0.1000e-03, 0.2000e-03, 0.3000e-03, 0.5000e-03,              
     & 0.7000e-03, 0.1100e-02, 0.1400e-02, 0.1800e-02, 0.2300e-02,              
     & 0.2700e-02, 0.3000e-02, 0.3300e-02, 0.3500e-02, 0.3700e-02,              
     & 0.3800e-02, 0.3800e-02, 0.3800e-02, 0.3600e-02, 0.3500e-02,              
     & 0.3300e-02, 0.3100e-02, 0.0000e+00, 0.1000e-03, 0.3000e-03,              
     & 0.4000e-03, 0.7000e-03, 0.9000e-03, 0.1300e-02, 0.1700e-02,              
     & 0.2000e-02, 0.2400e-02, 0.2700e-02, 0.3000e-02, 0.3200e-02,              
     & 0.3300e-02, 0.3400e-02, 0.3400e-02, 0.3300e-02, 0.3300e-02,              
     & 0.3100e-02, 0.3000e-02, 0.2800e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=13,15)/                                      
     & 0.0000e+00, 0.1000e-03, 0.2000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.8000e-03, 0.1100e-02, 0.1400e-02, 0.1800e-02, 0.2100e-02,              
     & 0.2400e-02, 0.2600e-02, 0.2700e-02, 0.2800e-02, 0.2900e-02,              
     & 0.3000e-02, 0.2900e-02, 0.2900e-02, 0.2800e-02, 0.2700e-02,              
     & 0.2500e-02, 0.0000e+00, 0.1000e-03, 0.1000e-03, 0.3000e-03,              
     & 0.4000e-03, 0.6000e-03, 0.9000e-03, 0.1100e-02, 0.1400e-02,              
     & 0.1700e-02, 0.2000e-02, 0.2200e-02, 0.2300e-02, 0.2400e-02,              
     & 0.2500e-02, 0.2600e-02, 0.2600e-02, 0.2600e-02, 0.2500e-02,              
     & 0.2400e-02, 0.2300e-02, 0.0000e+00, 0.0000e+00, 0.1000e-03,              
     & 0.2000e-03, 0.3000e-03, 0.4000e-03, 0.6000e-03, 0.9000e-03,              
     & 0.1100e-02, 0.1400e-02, 0.1600e-02, 0.1800e-02, 0.1900e-02,              
     & 0.2100e-02, 0.2200e-02, 0.2200e-02, 0.2300e-02, 0.2300e-02,              
     & 0.2300e-02, 0.2200e-02, 0.2200e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=16,18)/                                      
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.1000e-03, 0.2000e-03,              
     & 0.3000e-03, 0.4000e-03, 0.6000e-03, 0.8000e-03, 0.1000e-02,              
     & 0.1200e-02, 0.1400e-02, 0.1600e-02, 0.1700e-02, 0.1900e-02,              
     & 0.2000e-02, 0.2100e-02, 0.2100e-02, 0.2100e-02, 0.2100e-02,              
     & 0.2000e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.1000e-03, 0.2000e-03, 0.3000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.7000e-03, 0.9000e-03, 0.1100e-02, 0.1300e-02, 0.1500e-02,              
     & 0.1600e-02, 0.1800e-02, 0.1900e-02, 0.2000e-02, 0.2000e-02,              
     & 0.2000e-02, 0.1900e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.1000e-03, 0.1000e-03, 0.2000e-03,              
     & 0.4000e-03, 0.5000e-03, 0.7000e-03, 0.9000e-03, 0.1100e-02,              
     & 0.1200e-02, 0.1400e-02, 0.1600e-02, 0.1700e-02, 0.1800e-02,              
     & 0.1900e-02, 0.1900e-02, 0.1900e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=19,21)/                                      
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.1000e-03, 0.1000e-03, 0.2000e-03, 0.4000e-03,              
     & 0.5000e-03, 0.7000e-03, 0.9000e-03, 0.1100e-02, 0.1200e-02,              
     & 0.1400e-02, 0.1600e-02, 0.1700e-02, 0.1800e-02, 0.1800e-02,              
     & 0.1800e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.1000e-03, 0.1000e-03,              
     & 0.2000e-03, 0.4000e-03, 0.5000e-03, 0.7000e-03, 0.9000e-03,              
     & 0.1100e-02, 0.1300e-02, 0.1400e-02, 0.1600e-02, 0.1700e-02,              
     & 0.1700e-02, 0.1700e-02, 0.0000e+00,-0.1000e-03, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.1000e-03, 0.2000e-03, 0.3000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.8000e-03, 0.1000e-02, 0.1200e-02, 0.1300e-02, 0.1500e-02,              
     & 0.1600e-02, 0.1700e-02, 0.1700e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=22,24)/                                      
     &-0.1000e-03,-0.1000e-03, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.1000e-03,              
     & 0.2000e-03, 0.3000e-03, 0.5000e-03, 0.7000e-03, 0.9000e-03,              
     & 0.1100e-02, 0.1300e-02, 0.1400e-02, 0.1500e-02, 0.1600e-02,              
     & 0.1600e-02,-0.1000e-03,-0.1000e-03,-0.1000e-03, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.1000e-03, 0.1000e-03, 0.3000e-03, 0.4000e-03, 0.6000e-03,              
     & 0.8000e-03, 0.1000e-02, 0.1200e-02, 0.1400e-02, 0.1500e-02,              
     & 0.1600e-02, 0.1600e-02,-0.1000e-03,-0.1000e-03, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.1000e-03, 0.2000e-03, 0.4000e-03,              
     & 0.5000e-03, 0.7000e-03, 0.9000e-03, 0.1100e-02, 0.1300e-02,              
     & 0.1500e-02, 0.1600e-02, 0.1600e-02/                                      
      data ((o2(ip,iw),iw=1,21),ip=25,26)/                                      
     &-0.1000e-03,-0.1000e-03, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.1000e-03, 0.2000e-03, 0.3000e-03, 0.5000e-03, 0.7000e-03,              
     & 0.9000e-03, 0.1100e-02, 0.1300e-02, 0.1500e-02, 0.1600e-02,              
     & 0.1600e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,              
     & 0.0000e+00, 0.1000e-03, 0.2000e-03, 0.3000e-03, 0.5000e-03,              
     & 0.7000e-03, 0.9000e-03, 0.1100e-02, 0.1300e-02, 0.1500e-02,              
     & 0.1600e-02, 0.1600e-02/                                                  
c                                                                               
      data ((cah(i,j),i=1,22),j= 1, 5)/                                         
     & 0.9923, 0.9922, 0.9921, 0.9920, 0.9916, 0.9910, 0.9899, 0.9882,          
     & 0.9856, 0.9818, 0.9761, 0.9678, 0.9558, 0.9395, 0.9188, 0.8945,          
     & 0.8675, 0.8376, 0.8029, 0.7621, 0.7154, 0.6647, 0.9876, 0.9876,          
     & 0.9875, 0.9873, 0.9870, 0.9864, 0.9854, 0.9837, 0.9811, 0.9773,          
     & 0.9718, 0.9636, 0.9518, 0.9358, 0.9153, 0.8913, 0.8647, 0.8350,          
     & 0.8005, 0.7599, 0.7133, 0.6627, 0.9808, 0.9807, 0.9806, 0.9805,          
     & 0.9802, 0.9796, 0.9786, 0.9769, 0.9744, 0.9707, 0.9653, 0.9573,          
     & 0.9459, 0.9302, 0.9102, 0.8866, 0.8604, 0.8311, 0.7969, 0.7565,          
     & 0.7101, 0.6596, 0.9708, 0.9708, 0.9707, 0.9705, 0.9702, 0.9697,          
     & 0.9687, 0.9671, 0.9647, 0.9612, 0.9560, 0.9483, 0.9372, 0.9221,          
     & 0.9027, 0.8798, 0.8542, 0.8253, 0.7916, 0.7515, 0.7054, 0.6551,          
     & 0.9568, 0.9568, 0.9567, 0.9565, 0.9562, 0.9557, 0.9548, 0.9533,          
     & 0.9510, 0.9477, 0.9428, 0.9355, 0.9250, 0.9106, 0.8921, 0.8700,          
     & 0.8452, 0.8171, 0.7839, 0.7443, 0.6986, 0.6486/                          
      data ((cah(i,j),i=1,22),j= 6,10)/                                         
     & 0.9377, 0.9377, 0.9376, 0.9375, 0.9372, 0.9367, 0.9359, 0.9345,          
     & 0.9324, 0.9294, 0.9248, 0.9181, 0.9083, 0.8948, 0.8774, 0.8565,          
     & 0.8328, 0.8055, 0.7731, 0.7342, 0.6890, 0.6395, 0.9126, 0.9126,          
     & 0.9125, 0.9124, 0.9121, 0.9117, 0.9110, 0.9098, 0.9079, 0.9052,          
     & 0.9012, 0.8951, 0.8862, 0.8739, 0.8579, 0.8385, 0.8161, 0.7900,          
     & 0.7585, 0.7205, 0.6760, 0.6270, 0.8809, 0.8809, 0.8808, 0.8807,          
     & 0.8805, 0.8802, 0.8796, 0.8786, 0.8770, 0.8747, 0.8712, 0.8659,          
     & 0.8582, 0.8473, 0.8329, 0.8153, 0.7945, 0.7697, 0.7394, 0.7024,          
     & 0.6588, 0.6105, 0.8427, 0.8427, 0.8427, 0.8426, 0.8424, 0.8422,          
     & 0.8417, 0.8409, 0.8397, 0.8378, 0.8350, 0.8306, 0.8241, 0.8148,          
     & 0.8023, 0.7866, 0.7676, 0.7444, 0.7154, 0.6796, 0.6370, 0.5897,          
     & 0.7990, 0.7990, 0.7990, 0.7989, 0.7988, 0.7987, 0.7983, 0.7978,          
     & 0.7969, 0.7955, 0.7933, 0.7899, 0.7846, 0.7769, 0.7664, 0.7528,          
     & 0.7357, 0.7141, 0.6866, 0.6520, 0.6108, 0.5646/                          
      data ((cah(i,j),i=1,22),j=11,15)/                                         
     & 0.7515, 0.7515, 0.7515, 0.7515, 0.7514, 0.7513, 0.7511, 0.7507,          
     & 0.7501, 0.7491, 0.7476, 0.7450, 0.7409, 0.7347, 0.7261, 0.7144,          
     & 0.6992, 0.6793, 0.6533, 0.6203, 0.5805, 0.5357, 0.7020, 0.7020,          
     & 0.7020, 0.7019, 0.7019, 0.7018, 0.7017, 0.7015, 0.7011, 0.7005,          
     & 0.6993, 0.6974, 0.6943, 0.6894, 0.6823, 0.6723, 0.6588, 0.6406,          
     & 0.6161, 0.5847, 0.5466, 0.5034, 0.6518, 0.6518, 0.6518, 0.6518,          
     & 0.6518, 0.6517, 0.6517, 0.6515, 0.6513, 0.6508, 0.6500, 0.6485,          
     & 0.6459, 0.6419, 0.6359, 0.6273, 0.6151, 0.5983, 0.5755, 0.5458,          
     & 0.5095, 0.4681, 0.6017, 0.6017, 0.6017, 0.6017, 0.6016, 0.6016,          
     & 0.6016, 0.6015, 0.6013, 0.6009, 0.6002, 0.5989, 0.5967, 0.5932,          
     & 0.5879, 0.5801, 0.5691, 0.5535, 0.5322, 0.5043, 0.4700, 0.4308,          
     & 0.5518, 0.5518, 0.5518, 0.5518, 0.5518, 0.5518, 0.5517, 0.5516,          
     & 0.5514, 0.5511, 0.5505, 0.5493, 0.5473, 0.5441, 0.5393, 0.5322,          
     & 0.5220, 0.5076, 0.4878, 0.4617, 0.4297, 0.3929/                          
      data ((cah(i,j),i=1,22),j=16,19)/                                         
     & 0.5031, 0.5031, 0.5031, 0.5031, 0.5031, 0.5030, 0.5030, 0.5029,          
     & 0.5028, 0.5025, 0.5019, 0.5008, 0.4990, 0.4960, 0.4916, 0.4850,          
     & 0.4757, 0.4624, 0.4441, 0.4201, 0.3904, 0.3564, 0.4565, 0.4565,          
     & 0.4565, 0.4564, 0.4564, 0.4564, 0.4564, 0.4563, 0.4562, 0.4559,          
     & 0.4553, 0.4544, 0.4527, 0.4500, 0.4460, 0.4400, 0.4315, 0.4194,          
     & 0.4028, 0.3809, 0.3538, 0.3227, 0.4122, 0.4122, 0.4122, 0.4122,          
     & 0.4122, 0.4122, 0.4122, 0.4121, 0.4120, 0.4117, 0.4112, 0.4104,          
     & 0.4089, 0.4065, 0.4029, 0.3976, 0.3900, 0.3792, 0.3643, 0.3447,          
     & 0.3203, 0.2923, 0.3696, 0.3696, 0.3696, 0.3696, 0.3696, 0.3696,          
     & 0.3695, 0.3695, 0.3694, 0.3691, 0.3687, 0.3680, 0.3667, 0.3647,          
     & 0.3615, 0.3570, 0.3504, 0.3409, 0.3279, 0.3106, 0.2892, 0.2642/          
c                                                                               
      data (((lr(i,j,k),i=1,10),j=1,10),k= 1, 1)/                               
     &    338,    139,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    452,    236,    141,     91,     61,     42,          
     &     29,     20,     14,     10,    532,    359,    239,    164,          
     &    115,     82,     58,     42,     30,     22,    597,    477,          
     &    363,    273,    205,    154,    116,     87,     65,     48,          
     &    659,    569,    483,    401,    328,    265,    213,    169,          
     &    133,    105,    716,    645,    579,    516,    456,    398,          
     &    345,    297,    253,    215,    780,    725,    674,    627,          
     &    582,    538,    497,    457,    420,    384,    844,    805,          
     &    769,    736,    704,    674,    645,    618,    591,    566,          
     &    893,    866,    842,    819,    798,    777,    758,    740,          
     &    722,    705,    918,    897,    879,    862,    846,    831,          
     &    816,    803,    790,    778/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 2, 2)/                               
     &    338,    139,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    452,    236,    141,     91,     61,     42,          
     &     29,     20,     14,     10,    531,    359,    239,    164,          
     &    115,     82,     58,     42,     30,     22,    597,    477,          
     &    363,    273,    205,    154,    116,     86,     64,     48,          
     &    658,    569,    482,    400,    327,    265,    212,    169,          
     &    133,    104,    715,    644,    578,    515,    455,    398,          
     &    345,    296,    253,    214,    779,    723,    672,    625,          
     &    580,    536,    495,    455,    418,    382,    841,    802,          
     &    766,    732,    700,    670,    641,    613,    587,    562,          
     &    888,    861,    835,    812,    790,    769,    749,    731,          
     &    713,    696,    911,    889,    869,    851,    834,    817,          
     &    802,    788,    774,    761/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 3, 3)/                               
     &    337,    139,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    452,    236,    141,     91,     61,     42,          
     &     29,     20,     14,     10,    531,    358,    239,    164,          
     &    115,     82,     58,     42,     30,     22,    596,    476,          
     &    363,    273,    205,    154,    116,     86,     64,     48,          
     &    657,    568,    481,    399,    326,    264,    212,    168,          
     &    133,    104,    714,    642,    576,    513,    453,    396,          
     &    343,    295,    252,    214,    777,    721,    670,    622,          
     &    577,    534,    493,    453,    415,    380,    838,    798,          
     &    762,    728,    696,    665,    636,    608,    581,    556,          
     &    883,    854,    828,    804,    781,    759,    739,    720,          
     &    701,    684,    903,    879,    858,    838,    819,    802,          
     &    785,    770,    755,    741/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 4, 4)/                               
     &    337,    139,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    451,    236,    140,     91,     61,     42,          
     &     29,     20,     14,     10,    530,    358,    238,    164,          
     &    115,     82,     58,     42,     30,     22,    595,    476,          
     &    362,    272,    205,    154,    115,     86,     64,     48,          
     &    656,    567,    480,    398,    326,    263,    211,    168,          
     &    133,    104,    713,    641,    574,    512,    452,    395,          
     &    342,    294,    251,    213,    774,    718,    667,    619,          
     &    574,    531,    489,    450,    413,    377,    834,    794,          
     &    756,    722,    689,    659,    629,    601,    575,    549,          
     &    876,    846,    819,    794,    770,    748,    726,    707,          
     &    688,    670,    893,    867,    844,    823,    803,    784,          
     &    766,    750,    734,    719/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 5, 5)/                               
     &    337,    139,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    451,    236,    140,     91,     61,     42,          
     &     29,     20,     14,     10,    530,    358,    238,    164,          
     &    115,     82,     58,     42,     30,     22,    594,    475,          
     &    361,    272,    204,    154,    115,     86,     64,     48,          
     &    655,    565,    479,    397,    325,    263,    211,    168,          
     &    132,    104,    711,    639,    572,    510,    450,    393,          
     &    340,    293,    250,    212,    772,    715,    664,    616,          
     &    570,    527,    486,    446,    409,    374,    830,    788,          
     &    750,    715,    682,    651,    621,    593,    566,    541,          
     &    868,    836,    808,    781,    757,    733,    711,    691,          
     &    671,    653,    883,    854,    829,    806,    784,    764,          
     &    745,    727,    710,    694/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 6, 6)/                               
     &    337,    138,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    450,    235,    140,     90,     61,     41,          
     &     29,     20,     14,     10,    529,    357,    238,    163,          
     &    115,     82,     58,     42,     30,     22,    593,    474,          
     &    361,    271,    204,    153,    115,     86,     64,     48,          
     &    653,    564,    477,    396,    323,    262,    210,    167,          
     &    132,    103,    708,    636,    570,    507,    447,    391,          
     &    338,    291,    248,    210,    768,    711,    659,    611,          
     &    565,    522,    481,    442,    405,    370,    824,    781,          
     &    742,    706,    673,    641,    611,    583,    556,    531,          
     &    859,    825,    795,    767,    741,    717,    694,    672,          
     &    652,    633,    870,    839,    812,    786,    763,    741,          
     &    720,    701,    683,    666/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 7, 7)/                               
     &    336,    138,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    449,    235,    140,     90,     61,     41,          
     &     29,     20,     14,     10,    528,    356,    237,    163,          
     &    115,     81,     58,     42,     30,     22,    592,    473,          
     &    360,    270,    203,    153,    114,     86,     64,     48,          
     &    651,    561,    475,    394,    322,    260,    209,    166,          
     &    131,    103,    705,    633,    566,    504,    444,    388,          
     &    336,    288,    246,    209,    763,    706,    653,    605,          
     &    559,    516,    475,    436,    399,    365,    816,    772,          
     &    732,    696,    662,    630,    599,    571,    544,    518,          
     &    848,    812,    779,    750,    723,    697,    673,    651,          
     &    630,    610,    857,    823,    792,    765,    739,    716,          
     &    693,    673,    653,    636/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 8, 8)/                               
     &    336,    138,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    448,    234,    140,     90,     60,     41,          
     &     28,     20,     14,     10,    526,    355,    237,    163,          
     &    114,     81,     58,     41,     30,     22,    590,    471,          
     &    358,    269,    202,    152,    114,     85,     64,     48,          
     &    648,    559,    473,    392,    320,    259,    207,    165,          
     &    130,    102,    702,    629,    562,    500,    440,    384,          
     &    332,    285,    243,    206,    758,    699,    647,    598,          
     &    552,    509,    468,    429,    392,    358,    808,    762,          
     &    721,    683,    648,    615,    585,    556,    529,    503,          
     &    835,    796,    762,    730,    701,    674,    649,    626,          
     &    604,    584,    841,    804,    771,    741,    714,    688,          
     &    664,    642,    621,    603/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k= 9, 9)/                               
     &    335,    138,     76,     47,     31,     21,     14,     10,          
     &      7,      5,    447,    234,    139,     90,     60,     41,          
     &     28,     20,     14,     10,    524,    354,    236,    162,          
     &    114,     81,     58,     41,     30,     22,    587,    469,          
     &    357,    268,    201,    151,    113,     85,     63,     47,          
     &    645,    555,    470,    389,    318,    257,    206,    164,          
     &    129,    101,    697,    624,    557,    494,    435,    380,          
     &    328,    282,    240,    204,    751,    691,    638,    589,          
     &    543,    500,    459,    420,    384,    351,    797,    749,          
     &    707,    668,    632,    599,    567,    538,    511,    486,          
     &    820,    778,    741,    708,    677,    649,    622,    598,          
     &    575,    554,    825,    784,    748,    715,    686,    658,          
     &    632,    609,    587,    567/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=10,10)/                               
     &    334,    138,     76,     47,     31,     21,     14,      9,          
     &      7,      5,    445,    233,    139,     90,     60,     41,          
     &     28,     20,     14,     10,    522,    353,    235,    161,          
     &    113,     81,     57,     41,     30,     22,    584,    466,          
     &    354,    266,    200,    150,    113,     84,     63,     47,          
     &    641,    551,    466,    386,    315,    254,    204,    162,          
     &    128,    100,    691,    617,    551,    488,    429,    374,          
     &    323,    277,    236,    200,    743,    682,    628,    578,          
     &    532,    489,    448,    410,    374,    341,    784,    734,          
     &    690,    650,    613,    579,    547,    518,    490,    465,          
     &    803,    758,    718,    683,    650,    620,    592,    566,          
     &    543,    521,    806,    762,    723,    687,    655,    626,          
     &    599,    574,    551,    530/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=11,11)/                               
     &    333,    137,     76,     47,     31,     21,     14,      9,          
     &      7,      5,    443,    232,    138,     89,     60,     41,          
     &     28,     20,     14,     10,    519,    351,    234,    161,          
     &    113,     80,     57,     41,     29,     22,    581,    463,          
     &    352,    265,    199,    149,    112,     84,     62,     47,          
     &    636,    546,    461,    381,    311,    251,    201,    160,          
     &    126,     99,    684,    610,    543,    480,    422,    367,          
     &    317,    272,    232,    196,    732,    670,    615,    565,          
     &    519,    475,    435,    397,    362,    330,    770,    717,          
     &    671,    629,    591,    556,    524,    494,    466,    441,          
     &    784,    735,    693,    655,    620,    588,    559,    532,          
     &    508,    486,    786,    738,    695,    658,    623,    592,          
     &    563,    536,    513,    491/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=12,12)/                               
     &    331,    137,     76,     47,     31,     21,     14,      9,          
     &      7,      5,    441,    231,    138,     89,     60,     41,          
     &     28,     20,     14,     10,    516,    349,    232,    160,          
     &    112,     80,     57,     41,     29,     22,    576,    459,          
     &    349,    262,    197,    148,    111,     83,     62,     46,          
     &    629,    539,    455,    376,    307,    247,    198,    157,          
     &    124,     98,    676,    600,    533,    471,    413,    359,          
     &    310,    265,    226,    192,    720,    656,    600,    550,          
     &    503,    460,    420,    382,    348,    317,    752,    697,          
     &    648,    605,    566,    530,    497,    467,    439,    414,          
     &    763,    710,    665,    624,    587,    554,    523,    495,          
     &    470,    448,    764,    712,    666,    626,    589,    556,          
     &    525,    498,    473,    450/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=13,13)/                               
     &    329,    136,     75,     47,     31,     20,     14,      9,          
     &      7,      5,    438,    230,    137,     88,     59,     41,          
     &     28,     19,     14,     10,    512,    346,    230,    158,          
     &    111,     79,     56,     40,     29,     22,    570,    454,          
     &    345,    259,    194,    146,    109,     82,     61,     46,          
     &    622,    532,    447,    369,    301,    243,    194,    154,          
     &    122,     96,    665,    589,    521,    459,    402,    349,          
     &    300,    257,    219,    186,    705,    640,    582,    531,          
     &    484,    441,    401,    365,    331,    301,    732,    673,          
     &    623,    578,    537,    500,    466,    436,    408,    384,          
     &    740,    683,    634,    591,    552,    517,    485,    457,          
     &    431,    408,    740,    684,    635,    592,    553,    518,          
     &    486,    458,    432,    409/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=14,14)/                               
     &    327,    135,     75,     47,     30,     20,     14,      9,          
     &      7,      5,    434,    228,    136,     88,     59,     40,          
     &     28,     19,     14,     10,    506,    342,    228,    157,          
     &    110,     78,     56,     40,     29,     21,    563,    447,          
     &    340,    255,    192,    144,    108,     80,     60,     45,          
     &    612,    522,    438,    361,    294,    237,    190,    151,          
     &    119,     94,    652,    575,    507,    446,    389,    337,          
     &    290,    247,    211,    179,    688,    620,    562,    509,          
     &    462,    419,    380,    344,    312,    283,    710,    647,          
     &    594,    547,    505,    467,    433,    403,    375,    351,          
     &    715,    654,    602,    556,    515,    478,    446,    416,          
     &    390,    367,    715,    654,    602,    556,    515,    479,          
     &    446,    417,    391,    368/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=15,15)/                               
     &    324,    134,     74,     46,     30,     20,     14,      9,          
     &      6,      5,    429,    226,    135,     87,     58,     40,          
     &     27,     19,     14,     10,    499,    338,    225,    155,          
     &    109,     77,     55,     39,     28,     21,    554,    439,          
     &    334,    250,    188,    141,    106,     79,     59,     44,          
     &    600,    510,    427,    352,    286,    230,    184,    146,          
     &    115,     91,    637,    559,    491,    429,    373,    322,          
     &    277,    236,    201,    170,    668,    597,    537,    484,          
     &    437,    395,    356,    321,    290,    262,    684,    618,          
     &    562,    513,    470,    432,    397,    367,    340,    316,          
     &    687,    622,    567,    519,    476,    439,    405,    375,          
     &    349,    326,    687,    622,    567,    519,    477,    439,          
     &    405,    376,    349,    327/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=16,16)/                               
     &    320,    133,     74,     46,     30,     20,     13,      9,          
     &      6,      5,    423,    223,    133,     86,     58,     39,          
     &     27,     19,     13,     10,    491,    333,    222,    152,          
     &    107,     76,     54,     39,     28,     21,    543,    430,          
     &    326,    245,    183,    138,    103,     77,     58,     43,          
     &    586,    495,    414,    340,    276,    222,    177,    140,          
     &    111,     88,    618,    539,    471,    410,    355,    306,          
     &    262,    223,    189,    160,    644,    571,    510,    456,          
     &    409,    367,    329,    295,    266,    240,    656,    586,          
     &    528,    478,    433,    394,    360,    329,    303,    280,          
     &    658,    589,    531,    481,    437,    398,    364,    334,          
     &    308,    286,    658,    589,    531,    481,    437,    398,          
     &    364,    334,    308,    286/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=17,17)/                               
     &    316,    131,     73,     45,     30,     20,     13,      9,          
     &      6,      5,    416,    220,    131,     85,     57,     39,          
     &     27,     19,     13,     10,    481,    326,    217,    149,          
     &    105,     74,     53,     38,     27,     20,    529,    418,          
     &    317,    238,    178,    133,    100,     75,     56,     42,          
     &    568,    478,    398,    326,    264,    212,    169,    134,          
     &    106,     84,    597,    516,    448,    388,    334,    286,          
     &    244,    207,    176,    149,    617,    542,    479,    425,          
     &    378,    336,    299,    267,    239,    215,    625,    552,          
     &    491,    440,    394,    355,    321,    291,    265,    244,          
     &    626,    553,    493,    441,    396,    357,    323,    294,          
     &    268,    247,    626,    553,    493,    441,    396,    357,          
     &    323,    294,    268,    247/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=18,18)/                               
     &    310,    129,     72,     45,     29,     19,     13,      9,          
     &      6,      5,    407,    216,    129,     83,     56,     38,          
     &     26,     18,     13,     10,    468,    318,    212,    146,          
     &    102,     72,     52,     37,     27,     20,    513,    404,          
     &    306,    229,    172,    128,     96,     72,     54,     41,          
     &    548,    458,    379,    310,    250,    200,    159,    126,          
     &    100,     79,    571,    489,    421,    362,    310,    265,          
     &    225,    190,    161,    137,    587,    509,    445,    391,          
     &    344,    303,    268,    237,    211,    189,    592,    516,          
     &    453,    400,    355,    316,    282,    253,    229,    208,          
     &    593,    516,    454,    401,    356,    317,    283,    254,          
     &    230,    210,    593,    516,    454,    401,    356,    317,          
     &    283,    254,    230,    210/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=19,19)/                               
     &    303,    127,     71,     44,     29,     19,     13,      9,          
     &      6,      5,    396,    210,    126,     81,     55,     37,          
     &     26,     18,     13,      9,    453,    307,    205,    141,          
     &     99,     70,     50,     36,     26,     19,    494,    387,          
     &    293,    219,    164,    122,     91,     68,     51,     39,          
     &    523,    434,    357,    291,    234,    187,    148,    117,          
     &     93,     74,    542,    459,    392,    334,    284,    240,          
     &    203,    171,    145,    123,    553,    474,    409,    355,          
     &    309,    269,    235,    207,    183,    163,    556,    477,          
     &    413,    360,    315,    277,    244,    216,    193,    175,          
     &    557,    478,    414,    360,    315,    277,    244,    217,          
     &    194,    175,    557,    478,    414,    360,    315,    277,          
     &    244,    217,    194,    175/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=20,20)/                               
     &    295,    124,     69,     43,     28,     19,     13,      9,          
     &      6,      5,    382,    204,    122,     79,     53,     36,          
     &     25,     17,     12,      9,    434,    295,    197,    136,          
     &     95,     67,     48,     34,     25,     19,    470,    367,          
     &    277,    207,    154,    115,     86,     64,     48,     37,          
     &    495,    406,    332,    269,    215,    171,    136,    107,          
     &     85,     68,    509,    426,    359,    303,    255,    214,          
     &    180,    151,    127,    108,    517,    435,    371,    317,          
     &    272,    234,    203,    176,    154,    137,    518,    437,          
     &    373,    320,    276,    238,    207,    181,    160,    144,          
     &    519,    438,    373,    320,    276,    239,    207,    182,          
     &    161,    144,    519,    438,    373,    320,    276,    239,          
     &    207,    182,    161,    144/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=21,21)/                               
     &    284,    120,     67,     42,     27,     18,     12,      8,          
     &      6,      4,    365,    196,    118,     76,     51,     35,          
     &     24,     17,     12,      9,    412,    280,    188,    129,          
     &     90,     64,     46,     33,     24,     18,    443,    344,          
     &    259,    193,    144,    107,     80,     60,     45,     35,          
     &    462,    375,    304,    244,    195,    154,    122,     96,          
     &     76,     61,    473,    389,    324,    270,    225,    187,          
     &    156,    130,    109,     93,    477,    395,    331,    279,          
     &    236,    200,    171,    147,    127,    112,    478,    396,          
     &    332,    280,    237,    202,    173,    149,    130,    116,          
     &    478,    396,    332,    280,    237,    202,    173,    149,          
     &    130,    116,    478,    396,    332,    280,    237,    202,          
     &    173,    149,    130,    116/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=22,22)/                               
     &    271,    115,     64,     40,     26,     17,     12,      8,          
     &      6,      4,    345,    187,    112,     73,     49,     33,          
     &     23,     16,     11,      9,    386,    263,    176,    121,          
     &     85,     60,     43,     31,     22,     17,    411,    317,          
     &    237,    177,    131,     98,     73,     54,     41,     32,          
     &    426,    341,    273,    218,    173,    136,    107,     84,          
     &     67,     54,    433,    350,    287,    236,    194,    160,          
     &    132,    109,     91,     78,    435,    354,    291,    241,          
     &    200,    167,    140,    119,    102,     90,    436,    354,          
     &    291,    241,    201,    168,    141,    120,    104,     91,          
     &    436,    354,    291,    241,    201,    168,    141,    120,          
     &    104,     91,    436,    354,    291,    241,    201,    168,          
     &    141,    120,    104,     91/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=23,23)/                               
     &    255,    109,     61,     38,     25,     17,     11,      8,          
     &      5,      4,    321,    175,    105,     68,     46,     31,          
     &     21,     15,     11,      8,    355,    242,    163,    112,          
     &     78,     55,     39,     28,     21,     16,    375,    286,          
     &    214,    158,    117,     87,     65,     48,     37,     29,          
     &    385,    303,    240,    190,    149,    117,     91,     72,          
     &     57,     47,    390,    309,    249,    201,    163,    132,          
     &    108,     89,     74,     63,    391,    311,    251,    204,          
     &    166,    137,    113,     94,     80,     70,    391,    311,          
     &    251,    204,    167,    137,    113,     95,     81,     71,          
     &    391,    311,    251,    204,    167,    137,    113,     95,          
     &     81,     71,    391,    311,    251,    204,    167,    137,          
     &    113,     95,     81,     71/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=24,24)/                               
     &    236,    102,     57,     36,     23,     15,     10,      7,          
     &      5,      4,    293,    161,     97,     63,     42,     29,          
     &     20,     14,     10,      7,    320,    219,    147,    101,          
     &     70,     50,     35,     25,     19,     14,    335,    252,          
     &    188,    138,    102,     76,     56,     42,     32,     25,          
     &    341,    264,    206,    161,    125,     97,     75,     59,          
     &     47,     39,    344,    267,    211,    168,    134,    107,          
     &     86,     70,     58,     50,    344,    268,    212,    169,          
     &    135,    109,     88,     72,     61,     53,    344,    268,          
     &    212,    169,    135,    109,     88,     73,     61,     54,          
     &    344,    268,    212,    169,    135,    109,     88,     73,          
     &     61,     54,    344,    268,    212,    169,    135,    109,          
     &     88,     73,     61,     54/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=25,25)/                               
     &    214,     93,     52,     33,     21,     14,     10,      7,          
     &      5,      4,    260,    145,     88,     57,     38,     26,          
     &     18,     12,      9,      7,    281,    192,    129,     89,          
     &     62,     43,     31,     22,     16,     13,    291,    216,          
     &    159,    117,     86,     63,     47,     35,     27,     22,          
     &    295,    223,    171,    132,    101,     78,     60,     47,          
     &     37,     31,    296,    225,    174,    135,    106,     83,          
     &     66,     53,     44,     38,    296,    225,    174,    136,          
     &    106,     84,     67,     54,     45,     40,    296,    225,          
     &    174,    136,    106,     84,     67,     54,     45,     40,          
     &    296,    225,    174,    136,    106,     84,     67,     54,          
     &     45,     40,    296,    225,    174,    136,    106,     84,          
     &     67,     54,     45,     40/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=26,26)/                               
     &    187,     83,     47,     29,     19,     13,      9,      6,          
     &      4,      3,    224,    126,     77,     50,     33,     23,          
     &     16,     11,      8,      6,    238,    162,    109,     75,          
     &     52,     37,     26,     19,     14,     11,    244,    178,          
     &    130,     95,     69,     51,     37,     28,     22,     18,          
     &    246,    182,    137,    104,     79,     60,     45,     35,          
     &     28,     24,    246,    183,    138,    105,     81,     62,          
     &     48,     39,     32,     28,    246,    183,    138,    106,          
     &     81,     63,     49,     39,     32,     29,    246,    183,          
     &    138,    106,     81,     63,     49,     39,     32,     29,          
     &    246,    183,    138,    106,     81,     63,     49,     39,          
     &     32,     29,    246,    183,    138,    106,     81,     63,          
     &     49,     39,     32,     29/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=27,27)/                               
     &    156,     71,     40,     25,     16,     11,      7,      5,          
     &      4,      3,    183,    105,     64,     41,     28,     19,          
     &     13,      9,      7,      5,    192,    130,     88,     60,          
     &     42,     29,     21,     15,     11,      9,    195,    140,          
     &    101,     73,     53,     38,     28,     21,     16,     14,          
     &    195,    141,    104,     78,     58,     43,     33,     25,          
     &     20,     18,    196,    142,    105,     78,     59,     44,          
     &     34,     27,     22,     20,    196,    142,    105,     78,          
     &     59,     44,     34,     27,     22,     20,    196,    142,          
     &    105,     78,     59,     44,     34,     27,     22,     20,          
     &    196,    142,    105,     78,     59,     44,     34,     27,          
     &     22,     20,    196,    142,    105,     78,     59,     44,          
     &     34,     27,     22,     20/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=28,28)/                               
     &    122,     56,     32,     20,     13,      9,      6,      4,          
     &      3,      2,    139,     81,     50,     32,     21,     15,          
     &     10,      7,      5,      4,    143,     97,     65,     45,          
     &     31,     22,     15,     11,      8,      7,    145,    102,          
     &     73,     52,     37,     27,     19,     15,     11,     10,          
     &    145,    102,     74,     54,     39,     29,     22,     17,          
     &     14,     12,    145,    102,     74,     54,     40,     29,          
     &     22,     17,     14,     13,    145,    102,     74,     54,          
     &     40,     29,     22,     17,     14,     13,    145,    102,          
     &     74,     54,     40,     29,     22,     17,     14,     13,          
     &    145,    102,     74,     54,     40,     29,     22,     17,          
     &     14,     13,    145,    102,     74,     54,     40,     29,          
     &     22,     17,     14,     13/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=29,29)/                               
     &     83,     39,     23,     14,      9,      6,      4,      3,          
     &      2,      2,     93,     55,     34,     22,     15,     10,          
     &      7,      5,      3,      3,     94,     63,     43,     29,          
     &     20,     14,     10,      7,      5,      5,     95,     65,          
     &     46,     32,     23,     16,     12,      9,      7,      6,          
     &     95,     65,     46,     33,     24,     17,     13,     10,          
     &      8,      7,     95,     65,     46,     33,     24,     17,          
     &     13,     10,      8,      8,     95,     65,     46,     33,          
     &     24,     17,     13,     10,      8,      8,     95,     65,          
     &     46,     33,     24,     17,     13,     10,      8,      8,          
     &     95,     65,     46,     33,     24,     17,     13,     10,          
     &      8,      8,     95,     65,     46,     33,     24,     17,          
     &     13,     10,      8,      8/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=30,30)/                               
     &     42,     21,     12,      7,      5,      3,      2,      1,          
     &      1,      1,     46,     27,     17,     11,      7,      5,          
     &      3,      2,      2,      1,     46,     31,     21,     14,          
     &     10,      7,      5,      3,      3,      2,     46,     31,          
     &     22,     15,     10,      7,      5,      4,      3,      3,          
     &     46,     31,     22,     15,     11,      8,      5,      4,          
     &      3,      3,     46,     31,     22,     15,     11,      8,          
     &      5,      4,      4,      3,     46,     31,     22,     15,          
     &     11,      8,      5,      4,      4,      3,     46,     31,          
     &     22,     15,     11,      8,      5,      4,      4,      3,          
     &     46,     31,     22,     15,     11,      8,      5,      4,          
     &      4,      3,     46,     31,     22,     15,     11,      8,          
     &      5,      4,      4,      3/                                          
      data (((lr(i,j,k),i=1,10),j=1,10),k=31,31)/                               
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
c                                                                               
      data (((lt(i,j,k),i=1,10),j=1,10),k= 1, 1)/                               
     &    661,    860,    922,    952,    968,    979,    986,    990,          
     &    993,    995,    546,    762,    858,    908,    939,    958,          
     &    971,    980,    986,    989,    465,    639,    759,    834,          
     &    883,    917,    941,    957,    969,    977,    399,    519,          
     &    633,    724,    792,    843,    882,    912,    934,    950,          
     &    336,    425,    512,    594,    667,    730,    783,    827,          
     &    863,    893,    277,    347,    413,    475,    535,    593,          
     &    646,    695,    739,    778,    210,    263,    312,    358,          
     &    403,    445,    487,    526,    564,    600,    139,    175,          
     &    207,    238,    267,    296,    323,    349,    375,    399,          
     &     78,     98,    117,    134,    151,    167,    182,    197,          
     &    211,    225,     34,     43,     51,     58,     65,     72,          
     &     79,     85,     92,     98/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 2, 2)/                               
     &    660,    860,    922,    952,    968,    979,    986,    990,          
     &    993,    995,    545,    762,    858,    908,    938,    958,          
     &    971,    980,    986,    989,    465,    638,    759,    834,          
     &    883,    917,    940,    957,    969,    977,    398,    518,          
     &    633,    723,    792,    843,    882,    911,    934,    950,          
     &    336,    424,    511,    593,    667,    730,    783,    827,          
     &    863,    892,    276,    346,    412,    474,    534,    591,          
     &    645,    694,    738,    777,    209,    262,    310,    357,          
     &    401,    443,    484,    524,    562,    598,    137,    172,          
     &    205,    235,    264,    292,    319,    345,    371,    395,          
     &     76,     95,    113,    130,    146,    161,    176,    191,          
     &    205,    218,     31,     39,     46,     53,     60,     66,          
     &     72,     78,     84,     89/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 3, 3)/                               
     &    660,    860,    922,    952,    968,    979,    986,    990,          
     &    993,    995,    545,    762,    858,    908,    938,    958,          
     &    971,    980,    985,    989,    464,    638,    758,    834,          
     &    883,    917,    940,    957,    969,    977,    398,    518,          
     &    632,    723,    791,    842,    881,    911,    933,    950,          
     &    335,    423,    510,    592,    666,    729,    782,    826,          
     &    862,    891,    275,    345,    410,    472,    532,    590,          
     &    643,    692,    736,    776,    207,    260,    308,    354,          
     &    398,    441,    482,    521,    559,    595,    135,    170,          
     &    202,    232,    260,    288,    315,    341,    366,    390,          
     &     73,     92,    109,    125,    140,    155,    170,    184,          
     &    197,    210,     28,     35,     41,     47,     53,     59,          
     &     64,     70,     75,     80/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 4, 4)/                               
     &    660,    860,    922,    952,    968,    979,    986,    990,          
     &    993,    995,    544,    761,    858,    908,    938,    957,          
     &    971,    979,    985,    989,    464,    637,    758,    833,          
     &    882,    916,    940,    957,    969,    977,    397,    517,          
     &    631,    722,    790,    842,    881,    910,    933,    949,          
     &    334,    422,    509,    591,    664,    727,    780,    825,          
     &    861,    890,    274,    344,    408,    470,    530,    588,          
     &    641,    690,    734,    774,    205,    258,    306,    351,          
     &    395,    437,    478,    517,    555,    591,    133,    167,          
     &    198,    227,    256,    283,    309,    335,    359,    383,          
     &     69,     87,    103,    119,    134,    148,    162,    175,          
     &    188,    200,     24,     30,     36,     41,     47,     51,          
     &     56,     61,     65,     70/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 5, 5)/                               
     &    659,    859,    922,    951,    968,    979,    985,    990,          
     &    993,    995,    544,    761,    857,    908,    938,    957,          
     &    970,    979,    985,    989,    463,    636,    757,    833,          
     &    882,    916,    940,    956,    968,    976,    396,    516,          
     &    630,    721,    789,    841,    880,    910,    932,    949,          
     &    333,    421,    507,    589,    663,    726,    779,    823,          
     &    860,    889,    272,    342,    406,    468,    528,    585,          
     &    638,    687,    732,    771,    203,    255,    303,    348,          
     &    391,    433,    474,    513,    551,    587,    130,    163,          
     &    193,    222,    250,    277,    302,    328,    352,    376,          
     &     65,     82,     97,    112,    126,    139,    152,    165,          
     &    177,    189,     21,     26,     31,     35,     40,     44,          
     &     48,     52,     56,     59/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 6, 6)/                               
     &    658,    859,    922,    951,    968,    978,    985,    990,          
     &    993,    995,    543,    760,    857,    907,    937,    957,          
     &    970,    979,    985,    989,    462,    635,    756,    832,          
     &    881,    915,    939,    956,    968,    976,    395,    514,          
     &    629,    720,    788,    840,    879,    909,    931,    948,          
     &    331,    419,    505,    587,    661,    724,    777,    822,          
     &    858,    888,    270,    339,    403,    465,    525,    582,          
     &    635,    684,    729,    768,    200,    252,    299,    343,          
     &    386,    428,    468,    507,    545,    581,    126,    158,          
     &    188,    216,    243,    269,    294,    319,    343,    366,          
     &     61,     76,     90,    104,    117,    129,    141,    153,          
     &    165,    176,     17,     21,     25,     29,     33,     36,          
     &     39,     43,     46,     49/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 7, 7)/                               
     &    658,    859,    921,    951,    968,    978,    985,    990,          
     &    993,    994,    542,    760,    856,    907,    937,    957,          
     &    970,    979,    985,    989,    461,    634,    755,    831,          
     &    881,    915,    938,    955,    967,    976,    393,    512,          
     &    627,    718,    787,    839,    878,    908,    930,    947,          
     &    329,    417,    502,    585,    658,    721,    775,    819,          
     &    856,    886,    268,    336,    400,    461,    521,    578,          
     &    631,    680,    725,    765,    197,    247,    294,    338,          
     &    381,    422,    462,    501,    538,    574,    121,    152,          
     &    181,    208,    234,    260,    284,    308,    332,    354,          
     &     55,     69,     82,     95,    107,    118,    129,    140,          
     &    151,    161,     13,     17,     20,     23,     26,     29,          
     &     31,     34,     37,     39/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 8, 8)/                               
     &    656,    858,    921,    951,    968,    978,    985,    990,          
     &    993,    994,    540,    759,    855,    906,    937,    956,          
     &    969,    978,    984,    988,    459,    633,    754,    830,          
     &    880,    914,    938,    955,    967,    975,    391,    510,          
     &    625,    716,    785,    837,    876,    907,    929,    946,          
     &    327,    414,    499,    581,    655,    718,    772,    817,          
     &    854,    884,    265,    332,    396,    457,    516,    573,          
     &    626,    675,    720,    760,    193,    242,    288,    331,          
     &    373,    414,    454,    492,    529,    565,    116,    146,          
     &    173,    199,    224,    249,    272,    296,    318,    340,          
     &     49,     62,     74,     85,     95,    106,    116,    126,          
     &    135,    145,     10,     13,     15,     17,     20,     22,          
     &     24,     26,     28,     30/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k= 9, 9)/                               
     &    655,    857,    920,    950,    967,    978,    985,    989,          
     &    992,    994,    539,    757,    854,    905,    936,    956,          
     &    969,    978,    984,    988,    457,    631,    752,    828,          
     &    878,    913,    937,    954,    966,    974,    389,    507,          
     &    622,    713,    783,    835,    875,    905,    928,    945,          
     &    324,    410,    495,    577,    651,    715,    768,    813,          
     &    850,    881,    261,    328,    390,    451,    510,    566,          
     &    620,    669,    714,    754,    188,    236,    281,    323,          
     &    364,    404,    444,    482,    519,    555,    110,    138,          
     &    164,    188,    212,    235,    258,    281,    303,    324,          
     &     43,     54,     64,     74,     83,     93,    102,    110,          
     &    119,    127,      7,      9,     11,     13,     14,     16,          
     &     17,     19,     20,     22/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=10,10)/                               
     &    653,    856,    920,    950,    967,    977,    984,    989,          
     &    992,    994,    537,    756,    853,    904,    935,    955,          
     &    968,    977,    984,    987,    455,    628,    750,    827,          
     &    877,    911,    936,    953,    965,    974,    386,    504,          
     &    619,    710,    780,    832,    872,    903,    926,    943,          
     &    320,    406,    491,    572,    646,    710,    764,    809,          
     &    847,    877,    256,    322,    384,    444,    502,    559,          
     &    612,    661,    707,    747,    182,    229,    272,    313,          
     &    354,    393,    431,    469,    506,    542,    102,    128,          
     &    153,    176,    198,    220,    242,    263,    284,    305,          
     &     36,     46,     55,     63,     71,     79,     86,     94,          
     &    102,    109,      5,      6,      8,      9,     10,     11,          
     &     12,     13,     14,     15/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=11,11)/                               
     &    651,    855,    919,    949,    966,    977,    984,    989,          
     &    992,    994,    534,    754,    852,    903,    934,    954,          
     &    968,    977,    983,    987,    451,    625,    748,    825,          
     &    875,    910,    934,    952,    964,    972,    382,    500,          
     &    615,    706,    776,    829,    869,    900,    923,    941,          
     &    316,    400,    484,    566,    640,    704,    758,    804,          
     &    842,    873,    251,    315,    376,    435,    493,    549,          
     &    602,    652,    697,    738,    175,    220,    261,    301,          
     &    341,    379,    417,    454,    491,    526,     94,    118,          
     &    140,    162,    183,    203,    223,    243,    263,    283,          
     &     30,     38,     45,     52,     58,     65,     71,     78,          
     &     84,     90,      3,      4,      5,      6,      6,      7,          
     &      8,      8,      9,     10/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=12,12)/                               
     &    649,    854,    918,    948,    966,    977,    984,    988,          
     &    991,    993,    531,    751,    850,    902,    933,    953,          
     &    967,    976,    982,    986,    448,    621,    744,    822,          
     &    873,    908,    932,    950,    962,    971,    378,    494,          
     &    609,    701,    772,    825,    866,    897,    920,    938,          
     &    310,    394,    477,    558,    632,    697,    751,    797,          
     &    836,    867,    244,    307,    366,    424,    481,    537,          
     &    591,    640,    686,    727,    166,    209,    249,    287,          
     &    325,    362,    399,    436,    472,    507,     84,    106,          
     &    126,    146,    165,    184,    202,    221,    240,    258,          
     &     24,     30,     35,     41,     46,     52,     57,     62,          
     &     67,     72,      2,      2,      3,      3,      4,      4,          
     &      5,      5,      5,      6/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=13,13)/                               
     &    645,    852,    917,    947,    965,    976,    983,    988,          
     &    991,    993,    526,    748,    848,    900,    931,    952,          
     &    966,    975,    981,    985,    443,    616,    740,    818,          
     &    870,    905,    930,    948,    961,    969,    372,    488,          
     &    603,    695,    766,    820,    861,    892,    917,    935,          
     &    303,    385,    468,    549,    623,    688,    743,    789,          
     &    828,    860,    235,    296,    354,    411,    467,    523,          
     &    576,    626,    672,    714,    156,    196,    234,    271,          
     &    307,    343,    379,    414,    450,    485,     74,     93,          
     &    111,    128,    146,    163,    180,    197,    214,    231,          
     &     18,     22,     27,     31,     35,     39,     43,     47,          
     &     51,     56,      1,      1,      2,      2,      2,      2,          
     &      3,      3,      3,      3/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=14,14)/                               
     &    641,    850,    915,    946,    964,    975,    982,    987,          
     &    990,    992,    521,    744,    845,    898,    929,    950,          
     &    964,    974,    980,    984,    437,    610,    735,    814,          
     &    866,    902,    927,    945,    958,    967,    365,    480,          
     &    595,    688,    759,    813,    855,    887,    912,    931,          
     &    295,    375,    456,    537,    611,    676,    732,    779,          
     &    819,    852,    225,    284,    340,    395,    451,    506,          
     &    559,    609,    655,    698,    144,    182,    217,    252,          
     &    286,    320,    355,    389,    424,    459,     63,     80,          
     &     95,    110,    125,    140,    155,    171,    186,    202,          
     &     13,     16,     19,     22,     25,     28,     31,     34,          
     &     38,     41,      1,      1,      1,      1,      1,      1,          
     &      1,      1,      2,      2/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=15,15)/                               
     &    636,    847,    913,    945,    963,    974,    981,    986,          
     &    990,    992,    515,    739,    841,    895,    927,    948,          
     &    962,    972,    979,    983,    429,    603,    729,    808,          
     &    861,    897,    923,    942,    955,    965,    356,    470,          
     &    584,    678,    750,    806,    848,    881,    906,    925,          
     &    285,    362,    443,    523,    597,    662,    719,    767,          
     &    807,    841,    213,    269,    323,    376,    431,    485,          
     &    538,    588,    635,    678,    131,    165,    198,    230,          
     &    262,    294,    327,    361,    395,    428,     52,     66,          
     &     79,     92,    104,    117,    130,    144,    158,    172,          
     &      9,     11,     13,     15,     17,     19,     21,     24,          
     &     26,     28,      0,      0,      0,      0,      0,      1,          
     &      1,      1,      1,      1/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=16,16)/                               
     &    630,    843,    911,    943,    961,    973,    980,    985,          
     &    989,    991,    507,    733,    836,    891,    924,    945,          
     &    960,    970,    977,    981,    420,    594,    721,    802,          
     &    855,    892,    919,    938,    952,    961,    345,    457,          
     &    572,    667,    740,    796,    839,    873,    899,    918,          
     &    272,    348,    426,    505,    580,    646,    703,    752,          
     &    793,    828,    199,    252,    303,    354,    407,    460,          
     &    513,    563,    610,    654,    116,    147,    176,    206,          
     &    235,    265,    296,    328,    361,    394,     42,     52,          
     &     63,     73,     84,     95,    106,    117,    129,    141,          
     &      5,      7,      8,     10,     11,     12,     14,     15,          
     &     17,     18,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=17,17)/                               
     &    622,    839,    908,    940,    959,    971,    979,    984,          
     &    988,    990,    496,    726,    831,    886,    920,    942,          
     &    957,    967,    975,    979,    409,    582,    711,    793,          
     &    848,    886,    913,    933,    947,    957,    333,    442,          
     &    557,    652,    726,    784,    828,    863,    889,    910,          
     &    258,    330,    406,    484,    559,    625,    683,    733,          
     &    776,    812,    183,    232,    280,    329,    380,    432,          
     &    484,    534,    582,    626,    101,    127,    153,    179,          
     &    206,    234,    263,    293,    324,    356,     32,     40,          
     &     48,     56,     64,     73,     82,     92,    102,    112,          
     &      3,      4,      5,      6,      6,      7,      8,      9,          
     &     10,     11,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=18,18)/                               
     &    612,    833,    904,    937,    957,    969,    977,    982,          
     &    986,    989,    484,    716,    824,    881,    915,    938,          
     &    953,    964,    972,    977,    395,    568,    699,    783,          
     &    839,    878,    906,    926,    941,    952,    317,    424,          
     &    538,    635,    710,    769,    814,    850,    878,    900,          
     &    241,    309,    382,    460,    534,    601,    660,    711,          
     &    755,    792,    165,    210,    254,    300,    348,    399,          
     &    450,    500,    548,    593,     84,    107,    129,    152,          
     &    176,    200,    227,    255,    284,    315,     23,     29,          
     &     35,     41,     47,     54,     61,     68,     76,     85,          
     &      2,      2,      3,      3,      3,      4,      4,      5,          
     &      6,      6,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=19,19)/                               
     &    599,    826,    899,    934,    954,    966,    975,    980,          
     &    984,    987,    469,    704,    815,    873,    909,    933,          
     &    949,    960,    968,    973,    378,    551,    684,    770,          
     &    827,    868,    897,    918,    934,    945,    299,    402,          
     &    516,    613,    690,    751,    798,    835,    864,    887,          
     &    221,    285,    355,    431,    504,    572,    632,    684,          
     &    730,    769,    145,    185,    225,    268,    314,    362,          
     &    412,    461,    509,    555,     68,     87,    105,    124,          
     &    145,    166,    190,    215,    243,    271,     15,     19,          
     &     24,     28,     32,     37,     42,     48,     54,     61,          
     &      1,      1,      1,      1,      2,      2,      2,      2,          
     &      3,      3,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=20,20)/                               
     &    584,    817,    893,    929,    950,    963,    972,    978,          
     &    982,    985,    450,    690,    803,    864,    902,    926,          
     &    943,    955,    963,    969,    358,    529,    665,    753,          
     &    813,    855,    886,    908,    925,    937,    277,    376,          
     &    489,    588,    666,    728,    777,    816,    847,    871,          
     &    199,    257,    324,    397,    470,    538,    599,    652,          
     &    699,    740,    124,    159,    194,    233,    275,    322,          
     &    370,    418,    465,    511,     53,     67,     82,     97,          
     &    114,    133,    153,    176,    200,    226,     10,     12,          
     &     15,     18,     21,     24,     28,     32,     36,     42,          
     &      0,      0,      0,      1,      1,      1,      1,      1,          
     &      1,      1,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=21,21)/                               
     &    565,    806,    885,    923,    945,    959,    968,    975,          
     &    979,    982,    428,    672,    790,    853,    892,    918,          
     &    936,    949,    958,    964,    334,    504,    642,    734,          
     &    796,    840,    872,    896,    914,    927,    252,    346,          
     &    458,    557,    638,    702,    753,    793,    826,    852,          
     &    175,    227,    289,    359,    431,    499,    560,    615,          
     &    664,    706,    103,    132,    162,    196,    235,    278,          
     &    324,    370,    417,    462,     38,     49,     60,     73,          
     &     86,    101,    118,    137,    159,    182,      5,      7,          
     &      9,     10,     12,     14,     17,     19,     23,     26,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      1,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=22,22)/                               
     &    542,    793,    875,    915,    939,    954,    964,    971,          
     &    976,    979,    401,    650,    773,    839,    881,    908,          
     &    927,    941,    951,    958,    306,    474,    615,    710,          
     &    775,    821,    855,    881,    900,    915,    224,    312,          
     &    421,    521,    604,    670,    723,    766,    801,    829,          
     &    149,    194,    250,    317,    387,    454,    516,    572,          
     &    622,    666,     81,    105,    130,    159,    194,    233,          
     &    275,    320,    365,    409,     26,     34,     42,     51,          
     &     61,     73,     86,    102,    120,    141,      3,      4,          
     &      4,      5,      6,      8,      9,     11,     13,     15,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=23,23)/                               
     &    514,    776,    864,    906,    931,    947,    958,    966,          
     &    971,    975,    370,    623,    752,    822,    866,    895,          
     &    916,    931,    942,    950,    274,    438,    582,    681,          
     &    749,    798,    835,    862,    883,    899,    194,    273,          
     &    379,    480,    564,    632,    688,    733,    770,    800,          
     &    122,    161,    210,    272,    339,    405,    467,    523,          
     &    574,    620,     61,     79,    100,    124,    153,    187,          
     &    226,    267,    310,    352,     17,     21,     27,     33,          
     &     40,     49,     59,     72,     86,    103,      1,      2,          
     &      2,      2,      3,      4,      4,      5,      6,      8,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=24,24)/                               
     &    481,    755,    849,    895,    922,    939,    951,    959,          
     &    965,    970,    333,    591,    726,    802,    848,    880,          
     &    902,    919,    931,    940,    239,    396,    544,    646,          
     &    718,    771,    810,    839,    862,    880,    161,    232,          
     &    333,    433,    518,    589,    647,    694,    734,    766,          
     &     95,    127,    169,    225,    288,    352,    413,    469,          
     &    521,    567,     43,     57,     72,     91,    115,    144,          
     &    178,    215,    254,    294,      9,     12,     16,     20,          
     &     24,     30,     37,     46,     57,     70,      0,      1,          
     &      1,      1,      1,      1,      2,      2,      3,      4,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=25,25)/                               
     &    441,    729,    831,    881,    910,    929,    942,    951,          
     &    958,    964,    291,    552,    695,    776,    827,    861,          
     &    886,    904,    917,    927,    200,    349,    498,    606,          
     &    682,    737,    779,    812,    837,    857,    128,    188,          
     &    283,    380,    466,    539,    599,    649,    691,    726,          
     &     71,     95,    130,    178,    236,    296,    354,    410,          
     &    462,    509,     28,     37,     48,     62,     81,    104,          
     &    133,    165,    199,    235,      5,      6,      8,     10,          
     &     13,     17,     21,     27,     35,     45,      0,      0,          
     &      0,      0,      0,      0,      1,      1,      1,      1,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=26,26)/                               
     &    394,    698,    808,    863,    895,    916,    931,    942,          
     &    950,    956,    245,    507,    658,    745,    800,    838,          
     &    865,    885,    900,    912,    159,    297,    447,    558,          
     &    639,    698,    743,    779,    807,    829,     96,    145,          
     &    230,    324,    409,    483,    545,    597,    641,    679,          
     &     49,     67,     93,    133,    184,    239,    294,    348,          
     &    398,    445,     17,     22,     30,     39,     52,     70,          
     &     93,    119,    148,    180,      2,      3,      4,      5,          
     &      6,      8,     11,     14,     19,     26,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=27,27)/                               
     &    341,    661,    781,    841,    877,    901,    918,    930,          
     &    939,    946,    195,    454,    615,    708,    768,    810,          
     &    840,    863,    880,    893,    119,    241,    390,    504,          
     &    589,    652,    701,    739,    770,    795,     67,    103,          
     &    177,    265,    348,    422,    485,    539,    585,    625,          
     &     31,     43,     61,     93,    135,    184,    234,    284,          
     &    332,    377,      9,     12,     16,     22,     30,     43,          
     &     59,     80,    103,    129,      1,      1,      1,      2,          
     &      2,      3,      5,      7,      9,     13,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=28,28)/                               
     &    281,    616,    748,    815,    855,    882,    901,    915,          
     &    925,    933,    143,    395,    563,    665,    730,    776,          
     &    810,    835,    855,    870,     81,    184,    327,    444,          
     &    532,    599,    652,    693,    727,    755,     42,     67,          
     &    127,    206,    285,    357,    420,    475,    522,    564,          
     &     17,     24,     36,     59,     92,    132,    176,    221,          
     &    265,    308,      4,      5,      7,     10,     15,     23,          
     &     34,     49,     67,     87,      0,      0,      0,      1,          
     &      1,      1,      2,      3,      4,      6,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=29,29)/                               
     &    216,    564,    709,    783,    829,    859,    880,    896,          
     &    909,    918,     92,    329,    505,    614,    685,    736,          
     &    773,    802,    824,    842,     47,    128,    262,    378,          
     &    469,    539,    595,    640,    677,    708,     22,     37,          
     &     82,    150,    222,    290,    351,    406,    454,    496,          
     &      8,     11,     18,     32,     56,     88,    124,    163,          
     &    202,    240,      1,      2,      3,      4,      6,     10,          
     &     17,     27,     39,     54,      0,      0,      0,      0,          
     &      0,      0,      0,      1,      1,      2,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=30,30)/                               
     &    149,    503,    661,    745,    796,    831,    855,    874,          
     &    888,    900,     46,    260,    439,    555,    633,    689,          
     &    731,    763,    788,    809,     20,     78,    197,    309,          
     &    401,    474,    532,    580,    620,    653,      8,     15,          
     &     45,     99,    162,    224,    282,    335,    382,    424,          
     &      3,      4,      6,     14,     29,     52,     80,    111,          
     &    144,    177,      0,      0,      1,      1,      2,      4,          
     &      7,     13,     20,     30,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      1,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
      data (((lt(i,j,k),i=1,10),j=1,10),k=31,31)/                               
     &     82,    435,    607,    700,    757,    797,    825,    846,          
     &    863,    877,      7,    189,    368,    490,    574,    635,          
     &    681,    717,    745,    769,      0,     36,    135,    240,          
     &    329,    403,    463,    513,    555,    591,      0,      1,          
     &     18,     57,    108,    162,    215,    264,    308,    349,          
     &      0,      0,      0,      3,     12,     26,     46,     69,          
     &     95,    122,      0,      0,      0,      0,      0,      1,          
     &      2,      5,      9,     15,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0,      0,      0,      0,      0,          
     &      0,      0,      0,      0/                                          
c                                                                               
      data ((lr2(i,j),i=1,10),j=1,10)/                                          
     &    338,    139,     77,     48,     31,     21,     14,     10,          
     &      7,      5,    453,    237,    141,     91,     61,     42,          
     &     29,     20,     14,     10,    533,    360,    239,    165,          
     &    116,     82,     59,     42,     30,     22,    599,    479,          
     &    365,    274,    206,    155,    116,     87,     65,     48,          
     &    661,    572,    485,    403,    329,    267,    214,    170,          
     &    134,    105,    720,    649,    583,    520,    460,    402,          
     &    349,    300,    256,    217,    786,    732,    682,    634,          
     &    590,    546,    505,    465,    427,    391,    854,    817,          
     &    783,    751,    720,    691,    663,    635,    609,    584,          
     &    911,    888,    867,    848,    829,    811,    794,    777,          
     &    761,    746,    950,    937,    925,    914,    904,    894,          
     &    884,    875,    866,    857/                                          
      end                                                                       
  
********************* CLIRAD IR1  Date: October, 1994 ****************
*
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fito3 (npp1,p,at,ao)

      implicit none
c-----------------------------------------------------------------------
c--- Inputs: npp1 --- number of points in the vertical               ---
c---         p    --- pressure at the model levels                   ---
c-----------------------------------------------------------------------
c--- Outputs: at  --- tempreture at the model levels                 ---
c---          ao  --- ozone at the model levels                      ---
c-----------------------------------------------------------------------
      integer lay,n,npp1
      parameter (lay=88,n=61)

      real    p(lay),at(lay),ao(lay),pl(n),ta(n),wa(n),oa(n),pa(n)
      real    tl(n),wl(n),ol(n),w(n),tab(3),wk(n,4)
      integer itab(3),iop(2)
c----- ix=        1:trp; 2:mls; 3:mlw; 4:sas; 5:saw

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,ix,k,lun,int
      real    y

      save

      ix=1
      lun=90+ix

      OPEN( lun,FILE='trp.dat'
     1      ,FORM='FORMATTED',STATUS='OLD' )

      read (lun,51) (pl(i),tl(i),wl(i),ol(i),i=2,n)
   51 format (0p,2f10.3,1p,2e12.4)
      pl(1)=0.
      tl(1)=tl(2)
      wl(1)=wl(2)
      ol(1)=ol(2)
      do 25 i=1,n-1
        pa(i)=0.5*(pl(i)+pl(i+1))
        ta(i)=0.5*(tl(i)+tl(i+1))
        wa(i)=0.5*(wl(i)+wl(i+1))
        oa(i)=0.5*(ol(i)+ol(i+1))
   25 continue
      pa(n)=pl(n)
      ta(n)=tl(n)
      wa(n)=wl(n)
      oa(n)=ol(n)
      write (6,105)
      write (6,106) (i,pa(i),ta(i),oa(i),wa(i),i=1,n-1)

c-----------------------------------------------------------------------
c--- Starts interpolation                                            ---
c-----------------------------------------------------------------------
      iop(1)=4
      iop(2)=4
      int=1
      itab(1)=1
      itab(2)=0
      itab(3)=0
      call coeff (n,pa,ta,w,iop,int,wk)
      do 200 k=1,npp1-1
        y=p(k)
       call terp1 (n,pa,ta,w,y,int,tab,itab)
        at(k)=tab(1)
        print *,'k = ',k,'  y = ',y,'  at = ',at(k)
  200 continue
      call coeff (n,pa,oa,w,iop,int,wk)
      do 300 k=1,npp1-1
        y=p(k)
       call terp1 (n,pa,oa,w,y,int,tab,itab)
       ao(k)=tab(1)
       write(6,106) k,y,at(k),ao(k)
  300 continue
  105 format(5x,'k',7x,'pa',8x,'Ta',10x,'o3',8x,'qv')
  106 format(2x,i4,2f10.3,2e12.4)
      return
      end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine irrad (m,n,ndim,taucl,ccld,pl,ta,wa,oa,co2,ts,
     *                  high,flx,flc,dfdts,st4)

      implicit none

      integer NXX,nz,nadd,nbb,nxi,mm,np,n1,nx,no,nc,nh,nt
      parameter (NXX=514,NZ=43,nadd=7,nbb=1)
      parameter (nxi=nxx-2)
      parameter (mm=nxi/nbb/2,np=nz-2+nadd)
      parameter (n1=1)
      parameter (nx=26,no=21,nc=24,nh=31,nt=7)
c
c---- input parameters ------
      integer m,n,ndim
      real    co2
      real    taucl(m,ndim,np),ccld(m,ndim,np),pl(m,ndim,np+1),
     $        ta(m,ndim,np),wa(m,ndim,np),oa(m,ndim,np),ts(m,ndim)

      logical high

c---- output parameters ------

      real    flx(m,ndim,np+1),flc(m,ndim,np+1),dfdts(m,ndim,np+1),
     $        st4(m,ndim)

      real rflux(mm,7)
      common/radflux/rflux
*
**********************************************************************
*
* This routine computes ir fluxes due to water vapor, co2, and o3.
*   Clouds in different layers are assumed randomly overlapped.
*  
* This is a vectorized code.  It computes fluxes simultaneously for
*   (m x n) soundings, which is a subset of (m x ndim) soundings.
*   In a global climate model, m and ndim correspond to the numbers of
*   grid boxes in the zonal and meridional directions, respectively.
*
* Detailed description of the radiation routine is given in
*   Chou and Suarez (1994).
*
* There are two options for computing cooling rate profiles.
*
*   if high = .true., transmission functions in the co2, o3, and the
*   three water vapor bands with strong absorption are computed using
*   table look-up.  cooling rates are computed accurately from the
*   surface up to 0.01 mb.
*   if high = .false., transmission functions are computed using the
*   k-distribution method with linear pressure scaling.  cooling rates
*   are not calculated accurately for pressures less than 20 mb.
*   the computation is faster with high=.false. than with high=.true.
*
* The IR spectrum is divided into eight bands:
*   
*   bnad     wavenumber (/cm)   absorber         method
*
*    1           0 - 340           h2o            K/T
*    2         340 - 540           h2o            K/T
*    3         540 - 800       h2o,cont,co2       K,S,K/T
*    4         800 - 980       h2o,cont           K,S
*    5         980 - 1100      h2o,cont,o3        K,S,T
*    6        1100 - 1380      h2o,cont           K,S
*    7        1380 - 1900          h2o            K/T
*    8        1900 - 3000          h2o            K 
*
* Note : "h2o" for h2o line absorption
*        "cont" for h2o continuum absorption
*        "K" for k-distribution method
*        "S" for one-parameter temperature scaling
*        "T" for table look-up
*
* The 15 micrometer region (540-800/cm) is further divided into
*   3 sub-bands :
*
*   subbnad   wavenumber (/cm)
*
*    1          540 - 620
*    2          620 - 720
*    3          720 - 800
*
*---- Input parameters                               units    size
*
*   number of soundings in zonal direction (m)        n/d      1
*   number of soundings in meridional direction (n)   n/d      1
*   maximum number of soundings in
*                 meridional direction (ndim)         n/d      1 
*   number of atmospheric layers (np)                 n/d      1
*   cloud optical thickness (taucl)                   n/d     m*ndim*np
*   cloud cover (ccld)                              fraction  m*ndim*np
*   level pressure (pl)                               mb      m*ndim*(np+1)
*   layer temperature (ta)                            k       m*ndim*np
*   layer specific humidity (wa)                      g/g     m*ndim*np
*   layer ozone mixing ratio by mass (oa)             g/g     m*ndim*np
*   surface temperature (ts)                          k       m*ndim  
*   co2 mixing ratio by volumn (co2)                  pppv     1
*   high                                                       1
*
* pre-computed tables used in table look-up for transmittance calculations:
*
*   c1 , c2, c3: for co2 (band 3)
*   o1 , o2, o3: for  o3 (band 5)
*   h11,h12,h13: for h2o (band 1)
*   h21,h22,h23: for h2o (band 2)
*   h71,h72,h73: for h2o (band 7)
*
*---- output parameters
*
*   net downward flux, all-sky   (flx)             w/m**2     m*ndim*(np+1)
*   net downward flux, clear-sky (flc)             w/m**2     m*ndim*(np+1)
*   sensitivity of net downward flux  
*       to surface temperature (dfdts)             w/m**2/k   m*ndim*(np+1)
*   emission by the surface (st4)                  w/m**2     m*ndim 
* 
* Notes: 
*
*   (1)  Water vapor continuum absorption is included in 540-1380 /cm.
*   (2)  Scattering by clouds is not included.
*   (3)  Clouds are assumed "gray" bodies.
*   (4)  The diffuse cloud transmission is computed to be exp(-1.66*taucl).
*   (5)  If there are no clouds, flx=flc.
*   (6)  plevel(1) is the pressure at the top of the model atmosphere, and
*        plevel(np+1) is the surface pressure.
*   (7)  Downward flux is positive, and upward flux is negative.
*   (8)  dfdts is always negative because upward flux is defined as negative.
*   (9)  For questions and coding errors, please contact with Ming-Dah Chou,
*        Code 913, NASA/Goddard Space Flight Center, Greenbelt, MD 20771.
*        Phone: 301-286-4012, Fax: 301-286-1759,
*        e-mail: chou@climate.gsfc.nasa.gov
*
c-----parameters defining the size of the pre-computed tables for transmittance
c     calculations using table look-up.
c
c     "nx" is the number of intervals in pressure
c     "no" is the number of intervals in o3 amount
c     "nc" is the number of intervals in co2 amount
c     "nh" is the number of intervals in h2o amount
c     "nt" is the number of copies to be made from the pre-computed
c          transmittance tables to reduce "memory-bank conflict"
c          in parallel machines and, hence, enhancing the speed of
c          computations using table look-up. 
c          If such advantage does not exist, "nt" can be set to 1.
c***************************************************************************

ctao
      integer iradave
      common/iptionr/ iradave
c
c---- static data -----

      real    cb(5,8)

c---- temporary arrays -----

      real pa(MM,N1,np),dt(MM,N1,np)
      real sh2o(MM,N1,np+1),swpre(MM,N1,np+1),swtem(MM,N1,np+1)
      real sco3(MM,N1,np+1),scopre(MM,N1,np+1),scotem(MM,N1,np+1)
      real dh2o(MM,N1,np),dcont(MM,N1,np),dco2(MM,N1,np),do3(MM,N1,np)
      real th2o(MM,N1,6),tcon(MM,N1,3),tco2(MM,N1,6,2)
      real h2oexp(MM,N1,np,6), conexp(MM,N1,np,3),co2exp(MM,N1,np,6,2)
      real clr(MM,N1,0:np+1),fclr(MM,N1)
      real blayer(MM,N1,0:np+1),dbs(MM,N1)
      real trant(MM,N1)
      real flxu(MM,N1,np+1),flxd(MM,N1,np+1)

      logical oznbnd
      logical co2bnd
      logical h2otbl
      logical conbnd

      real c1 (nx,nc,nt),c2 (nx,nc,nt),c3 (nx,nc,nt)
      real o1 (nx,no,nt),o2 (nx,no,nt),o3 (nx,no,nt)
      real h11(nx,nh,nt),h12(nx,nh,nt),h13(nx,nh,nt)
      real h21(nx,nh,nt),h22(nx,nh,nt),h23(nx,nh,nt)
      real h71(nx,nh,nt),h72(nx,nh,nt),h73(nx,nh,nt)
      real dp,xx,w1,p1,dwe,dpe

c-----the following coefficients (table 2 of chou and suarez, 1995)
c     are for computing spectrally integtrated planck fluxes of
c     the 8 bands using eq. (22)
 
       data cb/
     1 -2.6844e-1,-8.8994e-2, 1.5676e-3,-2.9349e-6, 2.2233e-9,
     2  3.7315e+1,-7.4758e-1, 4.6151e-3,-6.3260e-6, 3.5647e-9,
     3  3.7187e+1,-3.9085e-1,-6.1072e-4, 1.4534e-5,-1.6863e-8,
     4 -4.1928e+1, 1.0027e+0,-8.5789e-3, 2.9199e-5,-2.5654e-8,
     5 -4.9163e+1, 9.8457e-1,-7.0968e-3, 2.0478e-5,-1.5514e-8,
     6 -1.0345e+2, 1.8636e+0,-1.1753e-2, 2.7864e-5,-1.1998e-8,
     7 -6.9233e+0,-1.5878e-1, 3.9160e-3,-2.4496e-5, 4.9301e-8,
     8  1.1483e+2,-2.2376e+0, 1.6394e-2,-5.3672e-5, 6.6456e-8/

c-----copy tables to enhance the speed of co2 (band 3), o3 (band5),
c     and h2o (bands 1, 2, and 7 only) transmission calculations
c     using table look-up.


      logical first
      data first /.true./
      integer i,j,k,ip,iw,it,ib,ik,iq,isb,k1,k2
      integer mmm


      include "h2o.tran3"
      include "co2.tran3"
      include "o3.tran3"

      save

ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1

      if (first) then

c-----tables co2 and h2o are only used with 'high' option

       if (high) then

        do iw=1,nh
         do ip=1,nx
          h11(ip,iw,1)=1.0-h11(ip,iw,1)
          h21(ip,iw,1)=1.0-h21(ip,iw,1)
          h71(ip,iw,1)=1.0-h71(ip,iw,1)
         enddo
        enddo

        do iw=1,nc
         do ip=1,nx
          c1(ip,iw,1)=1.0-c1(ip,iw,1)
         enddo
        enddo

c-----tables are replicated to avoid memory bank conflicts

        do it=2,nt
         do iw=1,nc
          do ip=1,nx
           c1 (ip,iw,it)= c1(ip,iw,1)
           c2 (ip,iw,it)= c2(ip,iw,1)
           c3 (ip,iw,it)= c3(ip,iw,1)
          enddo
         enddo
         do iw=1,nh
          do ip=1,nx
           h11(ip,iw,it)=h11(ip,iw,1)
           h12(ip,iw,it)=h12(ip,iw,1)
           h13(ip,iw,it)=h13(ip,iw,1)
           h21(ip,iw,it)=h21(ip,iw,1)
           h22(ip,iw,it)=h22(ip,iw,1)
           h23(ip,iw,it)=h23(ip,iw,1)
           h71(ip,iw,it)=h71(ip,iw,1)
           h72(ip,iw,it)=h72(ip,iw,1)
           h73(ip,iw,it)=h73(ip,iw,1)
          enddo
         enddo
        enddo

       endif

c-----always use table look-up for ozone transmittance

        do iw=1,no      
         do ip=1,nx
          o1(ip,iw,1)=1.0-o1(ip,iw,1)
         enddo
        enddo

        do it=2,nt
         do iw=1,no
          do ip=1,nx
           o1 (ip,iw,it)= o1(ip,iw,1)
           o2 (ip,iw,it)= o2(ip,iw,1)
           o3 (ip,iw,it)= o3(ip,iw,1)
          enddo
         enddo
        enddo

       first=.false.

      endif

c-----compute layer pressure (pa) and layer temperature minus 250K (dt)
 
      do k=1,np
       do j=1,n
        do i=1,mmm
         pa(i,j,k)=0.5*(pl(i,j,k)+pl(i,j,k+1))
         dt(i,j,k)=ta(i,j,k)-250.0
        enddo
       enddo
      enddo

c-----compute layer absorber amount

c     dh2o : water vapor amount (g/cm**2)
c     dcont: scaled water vapor amount for continuum absorption (g/cm**2)
c     dco2 : co2 amount (cm-atm)stp
c     do3  : o3 amount (cm-atm)stp
c     the factor 1.02 is equal to 1000/980
c     factors 789 and 476 are for unit conversion
c     the factor 0.001618 is equal to 1.02/(.622*1013.25) 
c     the factor 6.081 is equal to 1800/296

      do k=1,np
       do j=1,n
        do i=1,mmm

         dp           = pl(i,j,k+1)-pl(i,j,k)
         dh2o (i,j,k) = 1.02*wa(i,j,k)*dp+1.e-10
         dco2 (i,j,k) = 789.*co2*dp+1.e-10
         do3  (i,j,k) = 476.0*oa(i,j,k)*dp+1.e-10

c-----compute scaled water vapor amount for h2o continuum absorption
c     following eq. (43).

         xx=pa(i,j,k)*0.001618*wa(i,j,k)*wa(i,j,k)*dp
         dcont(i,j,k) = xx*exp(1800./ta(i,j,k)-6.081)+1.e-10

c-----compute effective cloud-free fraction, clr, for each layer.
c     the cloud diffuse transmittance is approximated by using a
c     diffusivity factor of 1.66.

         clr(i,j,k)=1.0-(ccld(i,j,k)*(1.-exp(-1.66*taucl(i,j,k))))

        enddo
       enddo
      enddo

c-----compute column-integrated h2o amoumt, h2o-weighted pressure
c     and temperature.  it follows eqs. (37) and (38).

       if (high) then

        call column(m,n,np,pa,dt,dh2o,sh2o,swpre,swtem)

       endif

c-----the surface (with an index np+1) is treated as a layer filled with
c     black clouds.

      do j=1,n
       do i=1,mmm
        clr(i,j,0)    = 1.0
        clr(i,j,np+1) = 0.0
        st4(i,j)      = 0.0
       enddo
      enddo

c-----initialize fluxes

      do k=1,np+1
       do j=1,n
        do i=1,mmm
         flx(i,j,k)  = 0.0
         flc(i,j,k)  = 0.0
         dfdts(i,j,k)= 0.0
         flxu(i,j,k) = 0.0
         flxd(i,j,k) = 0.0
        enddo
       enddo
      enddo

c-----integration over spectral bands

      do 1000 ib=1,8

c-----if h2otbl, compute h2o (line) transmittance using table look-up.
c     if conbnd, compute h2o (continuum) transmittance in bands 3, 4, 5 and 6.
c     if co2bnd, compute co2 transmittance in band 3.
c     if oznbnd, compute  o3 transmittance in band 5.
       h2otbl=high.and.(ib.eq.1.or.ib.eq.2.or.ib.eq.7)
       conbnd=ib.ge.3.and.ib.le.6
       co2bnd=ib.eq.3
       oznbnd=ib.eq.5

c-----blayer is the spectrally integrated planck flux of the mean layer
c     temperature derived from eq. (22)
c     the fitting for the planck flux is valid in the range 160-345 K.

       do k=1,np
        do j=1,n
         do i=1,mmm
          blayer(i,j,k)=ta(i,j,k)*(ta(i,j,k)*(ta(i,j,k)
     *                 *(ta(i,j,k)*cb(5,ib)+cb(4,ib))+cb(3,ib))
     *                 +cb(2,ib))+cb(1,ib)
         enddo
        enddo
       enddo

c-----the earth's surface, with an index "np+1", is treated as a layer

       do j=1,n
        do i=1,mmm
         blayer(i,j,0)   = 0.0
         blayer(i,j,np+1)=ts(i,j)*(ts(i,j)*(ts(i,j)
     *                   *(ts(i,j)*cb(5,ib)+cb(4,ib))+cb(3,ib))
     *                   +cb(2,ib))+cb(1,ib)

c-----dbs is the derivative of the surface planck flux with respect to
c     surface temperature (eq. 59).

         dbs(i,j)=ts(i,j)*(ts(i,j)*(ts(i,j)
     *           *4.*cb(5,ib)+3.*cb(4,ib))+2.*cb(3,ib))+cb(2,ib)

        enddo
       enddo

c-----compute column-integrated absorber amoumt, absorber-weighted
c     pressure and temperature for co2 (band 3) and o3 (band 5).
c     it follows eqs. (37) and (38).

c-----this is in the band loop to save storage

      if( high .and. co2bnd) then

        call column(m,n,np,pa,dt,dco2,sco3,scopre,scotem)

      endif

      if(oznbnd) then

        call column(m,n,np,pa,dt,do3,sco3,scopre,scotem)

      endif

c-----compute the exponential terms (eq. 32) at each layer for
c     water vapor line absorption when k-distribution is used

      if( .not. h2otbl) then

        call h2oexps(ib,m,n,np,dh2o,pa,dt,h2oexp)

      endif

c-----compute the exponential terms (eq. 46) at each layer for
c     water vapor continuum absorption

      if( conbnd) then

        call conexps(ib,m,n,np,dcont,conexp)

      endif


c-----compute the  exponential terms (eq. 32) at each layer for
c     co2 absorption

      if( .not.high .and. co2bnd) then

        call co2exps(m,n,np,dco2,pa,dt,co2exp)

      endif

c-----compute transmittances for regions between levels k1 and k2
c     and update the fluxes at the two levels.

      do 2000 k1=1,np

c-----initialize fclr, th2o, tcon, and tco2

        do j=1,n
         do i=1,mmm
          fclr(i,j)=1.0
         enddo
        enddo

c-----for h2o line absorption

      if(.not. h2otbl) then
        do ik=1,6
         do j=1,n
          do i=1,mmm
           th2o(i,j,ik)=1.0
          enddo
         enddo
        enddo
      endif

c-----for h2o continuum absorption

       if (conbnd) then
         do iq=1,3
          do j=1,n
           do i=1,mmm
            tcon(i,j,iq)=1.0
           enddo
          enddo
         enddo
       endif

c-----for co2 absorption when using k-distribution method.
c     band 3 is divided into 3 sub-bands, but sub-bands 3a and 3c
c     are combined in computing the co2 transmittance.

       if (.not. high .and. co2bnd) then
         do isb=1,2
          do ik=1,6
           do j=1,n
            do i=1,mmm
             tco2(i,j,ik,isb)=1.0
            enddo
           enddo
          enddo
         enddo
       endif

c-----loop over the bottom level of the region (k2)

      do 3000 k2=k1+1,np+1

          do j=1,n
           do i=1,mmm
            trant(i,j)=1.0
           enddo
          enddo

       if(h2otbl) then

          w1=-8.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

c-----compute water vapor transmittance using table look-up

          if (ib.eq.1 ) then

           call tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,
     *                 w1,p1,dwe,dpe,h11,h12,h13,trant)

          endif
          if (ib.eq.2 ) then

           call tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,
     *                 w1,p1,dwe,dpe,h21,h22,h23,trant)

          endif
          if (ib.eq.7 ) then

           call tablup(k1,k2,m,n,np,nx,nh,nt,sh2o,swpre,swtem,
     *                 w1,p1,dwe,dpe,h71,h72,h73,trant)

          endif

       else

c-----compute water vapor transmittance using k-distribution.

        call wvkdis(ib,m,n,np,k2-1,h2oexp,conexp,th2o,tcon,trant)

       endif

       if(co2bnd) then

        if( high ) then

c-----compute co2 transmittance using table look-up method

          w1=-4.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

          call tablup(k1,k2,m,n,np,nx,nc,nt,sco3,scopre,scotem,
     *                w1,p1,dwe,dpe,c1,c2,c3,trant)

        else

c-----compute co2 transmittance using k-distribution method

          call co2kdis(m,n,np,k2-1,co2exp,tco2,trant)
        
        endif 

       endif 

c-----compute o3 transmittance using table look-up

       if (oznbnd) then

          w1=-6.0
          p1=-2.0
          dwe=0.3
          dpe=0.2

          call tablup(k1,k2,m,n,np,nx,no,nt,sco3,scopre,scotem,
     *                w1,p1,dwe,dpe,o1,o2,o3,trant)

       endif

c-----fclr is the clear line-of-sight between levels k1 and k2.
c     in computing fclr, clouds are assumed randomly overlapped
c     using eq. (10).
 
      do j=1,n
       do i=1,mmm
        fclr(i,j) = fclr(i,j)*clr(i,j,k2-1)
       enddo
      enddo

c-----compute upward and downward fluxes


c-----add "boundary" terms to the net downward flux.
c     these are the first terms on the right-hand-side of
c     eqs. (56a) and (56b).
c     downward fluxes are positive.

      if (k2 .eq. k1+1) then
       do j=1,n
        do i=1,mmm
         flc(i,j,k1)=flc(i,j,k1)-blayer(i,j,k1)
         flc(i,j,k2)=flc(i,j,k2)+blayer(i,j,k1)
        enddo
       enddo
      endif

c-----add flux components involving the four layers above and below
c     the levels k1 and k2.  it follows eqs. (56a) and (56b).

      do j=1,n
       do i=1,mmm
        xx=trant(i,j)*(blayer(i,j,k2-1)-blayer(i,j,k2))
        flc(i,j,k1) =flc(i,j,k1)+xx
        xx=trant(i,j)*(blayer(i,j,k1-1)-blayer(i,j,k1))
        flc(i,j,k2) =flc(i,j,k2)+xx
       enddo
      enddo

c-----compute upward and downward fluxes for all-sky situation

      if (k2 .eq. k1+1) then
       do j=1,n
        do i=1,mmm
         flxu(i,j,k1)=flxu(i,j,k1)-blayer(i,j,k1)
         flxd(i,j,k2)=flxd(i,j,k2)+blayer(i,j,k1)
        enddo
       enddo
      endif

      do j=1,n
       do i=1,mmm
        xx=trant(i,j)*(blayer(i,j,k2-1)-blayer(i,j,k2))
        flxu(i,j,k1) =flxu(i,j,k1)+xx*fclr(i,j)
        xx=trant(i,j)*(blayer(i,j,k1-1)-blayer(i,j,k1))
        flxd(i,j,k2) =flxd(i,j,k2)+xx*fclr(i,j)
       enddo
      enddo


 3000 continue

c-----compute the partial derivative of fluxes with respect to
c     surface temperature (eq. 59).

      do j=1,n
       do i=1,mmm
        dfdts(i,j,k1) =dfdts(i,j,k1)-dbs(i,j)*trant(i,j)*fclr(i,j)
       enddo
      enddo
     
 2000 continue

c-----add contribution from the surface to the flux terms at the surface.

      do j=1,n
       do i=1,mmm
        dfdts(i,j,np+1) =dfdts(i,j,np+1)-dbs(i,j)
       enddo
      enddo

      do j=1,n
       do i=1,mmm
        flc(i,j,np+1)=flc(i,j,np+1)-blayer(i,j,np+1)
        flxu(i,j,np+1)=flxu(i,j,np+1)-blayer(i,j,np+1)
        st4(i,j)=st4(i,j)-blayer(i,j,np+1)
       enddo
      enddo


c     write(7,3211) ib, flxd(1,1,52),flxu(1,1,52)
c     write(7,3211) ib, flxd(1,1,np+1),flxu(1,1,np+1)
c     3211 format ('ib, fluxd, fluxu=', i3,2f12.3)

 1000 continue

      do k=1,np+1
       do j=1,n
        do i=1,mmm
         flx(i,j,k)=flxd(i,j,k)+flxu(i,j,k)
        enddo
       enddo
      enddo


C --------------------  GCSS Wrokshop -----------------------------------

      do i=1,mm
        rflux(i,2)=flxd(i,1,np+1)
        rflux(i,5)=flxu(i,1,1)
        rflux(i,7)=flxu(i,1,np+1)
      enddo

C -----------------------------------------------------------------------

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine column (m,n,np,pa,dt,sabs0,sabs,spre,stem)
c
c**************************************************************************
c-----compute column-integrated (from top of the model atmosphere)
c     absorber amount (sabs), absorber-weighted pressure (spre) and
c     temperature (stem).
c     computations of spre and stem follows eqs. (37) and (38).
c
c--- input parameters
c   number of soundings in zonal direction (m)
c   number of soundings in meridional direction (n)
c   number of atmospheric layers (np)
c   layer pressure (pa)
c   layer temperature minus 250K (dt)
c   layer absorber amount (sabs0)
c
c--- output parameters
c   column-integrated absorber amount (sabs)
c   column absorber-weighted pressure (spre)
c   column absorber-weighted temperature (stem)
c
c--- units of pa and dt are mb and k, respectively.
c    units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3
c**************************************************************************
      implicit none
      integer m,n,np

c---- input parameters -----

      real    pa(m,n,np),dt(m,n,np),sabs0(m,n,np)

c---- output parameters -----

      real    sabs(m,n,np+1),spre(m,n,np+1),stem(m,n,np+1)

      integer iradave
      common/iptionr/ iradave

      integer i,j,k,mmm

      save

c*********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
        do j=1,n
         do i=1,mmm
          sabs(i,j,1)=0.0
          spre(i,j,1)=0.0
          stem(i,j,1)=0.0
         enddo
        enddo

        do k=1,np
         do j=1,n
          do i=1,mmm
           sabs(i,j,k+1)=sabs(i,j,k)+sabs0(i,j,k)
           spre(i,j,k+1)=spre(i,j,k)+pa(i,j,k)*sabs0(i,j,k)
           stem(i,j,k+1)=stem(i,j,k)+dt(i,j,k)*sabs0(i,j,k)
          enddo
         enddo
        enddo

       return
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tablup(k1,k2,m,n,np,nx,nh,nt,sabs,spre,stem,w1,p1,
     *                  dwe,dpe,coef1,coef2,coef3,tran)
c**********************************************************************
c   compute water vapor, co2, and o3 transmittances between levels k1 and k2
c   using table look-up for m x n soundings.
c
c   Calculations follow Eq. (40) of Chou and Suarez (1995)
c
c---- input ---------------------
c  indices for pressure levels (k1 and k2)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of atmospheric layers (np)
c  number of pressure intervals in the table (nx)
c  number of absorber amount intervals in the table (nh)
c  number of tables copied (nt)
c  column-integrated absorber amount (sabs)
c  column absorber amount-weighted pressure (spre)
c  column absorber amount-weighted temperature (stem)
c  first value of absorber amount (log10) in the table (w1) 
c  first value of pressure (log10) in the table (p1) 
c  size of the interval of absorber amount (log10) in the table (dwe)
c  size of the interval of pressure (log10) in the table (dpe)
c  pre-computed coefficients (coef1, coef2, and coef3)
c
c---- updated ---------------------
c  transmittance (tran)
c
c  Note:
c   (1) units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3.
c   (2) units of spre and stem are, respectively, mb and K.
c   (3) there are nt identical copies of the tables (coef1, coef2, and
c       coef3).  the prupose of using the multiple copies of tables is
c       to increase the speed in parallel (vectorized) computations.
C       if such advantage does not exist, nt can be set to 1.
c   
c**********************************************************************
      implicit none
      integer k1,k2,m,n,np,nx,nh,nt

c---- input parameters -----

      real    w1,p1,dwe,dpe
      real    sabs(m,n,np+1),spre(m,n,np+1),stem(m,n,np+1)
      real    coef1(nx,nh,nt),coef2(nx,nh,nt),coef3(nx,nh,nt)

c---- update parameter -----

      real tran(m,n)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- temporary variables -----

      integer i,j,mmm,iw,ip,nn
      real    x1,x2,x3,we,pe,fw,fp,pa,pb,pc,ax,ba,bb,t1,ca,cb,t2

      save

c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
      do j=1,n
       do i=1,mmm

        nn=mod(i,nt)+1

        x1=sabs(i,j,k2)-sabs(i,j,k1)
        x2=(spre(i,j,k2)-spre(i,j,k1))/x1
        x3=(stem(i,j,k2)-stem(i,j,k1))/x1

        we=(log10(x1)-w1)/dwe
        pe=(log10(x2)-p1)/dpe

        we=max(we,w1-2.*dwe)
        pe=max(pe,p1)

        iw=int(we+1.5)
        ip=int(pe+1.5)

        iw=min(iw,nh-1)
        iw=max(iw, 2)

        ip=min(ip,nx-1)
        ip=max(ip, 1)

        fw=we-float(iw-1)
        fp=pe-float(ip-1)

c-----linear interpolation in pressure

        pa = coef1(ip,iw-1,nn)*(1.-fp)+coef1(ip+1,iw-1,nn)*fp
        pb = coef1(ip,iw,  nn)*(1.-fp)+coef1(ip+1,iw,  nn)*fp
        pc = coef1(ip,iw+1,nn)*(1.-fp)+coef1(ip+1,iw+1,nn)*fp

c-----quadratic interpolation in absorber amount for coef1

        ax = (-pa*(1.-fw)+pc*(1.+fw)) *fw*0.5 + pb*(1.-fw*fw)

c-----linear interpolation in absorber amount for coef2 and coef3

        ba = coef2(ip,iw,  nn)*(1.-fp)+coef2(ip+1,iw,  nn)*fp
        bb = coef2(ip,iw+1,nn)*(1.-fp)+coef2(ip+1,iw+1,nn)*fp
        t1 = ba*(1.-fw) + bb*fw

        ca = coef3(ip,iw,  nn)*(1.-fp)+coef3(ip+1,iw,  nn)*fp
        cb = coef3(ip,iw+1,nn)*(1.-fp)+coef3(ip+1,iw+1,nn)*fp
        t2 = ca*(1.-fw) + cb*fw

c-----update the total transmittance between levels k1 and k2

        tran(i,j)= (ax + (t1+t2*x3) * x3)*tran(i,j)

       enddo
      enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wvkdis(ib,m,n,np,k,h2oexp,conexp,th2o,tcon,tran)
c**********************************************************************
c   compute water vapor transmittance between levels k1 and k2 for
c   m x n soundings using the k-distribution method.
c
c   computations follow eqs. (34), (46), (50) and (52).
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of levels (np)
c  current level (k)
c  exponentials for line absorption (h2oexp) 
c  exponentials for continuum absorption (conexp) 
c
c---- updated parameters
c  transmittance between levels k1 and k2 due to
c    water vapor line absorption (th2o)
c  transmittance between levels k1 and k2 due to
c    water vapor continuum absorption (tcon)
c  total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np,k

c---- input parameters ------

      real    conexp(m,n,np,3),h2oexp(m,n,np,6)

c---- updated parameters -----

      real    th2o(m,n,6),tcon(m,n,3),tran(m,n)

      integer iradave
      common/iptionr/ iradAVE

c---- static data -----
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------

      integer ne(8)
      real    fkw(6,8),gkw(6,3)

c---- temporary arrays -----
      real    trnth2o

c-----fkw is the planck-weighted k-distribution function due to h2o
c     line absorption given in table 4 of Chou and Suarez (1995).
c     the k-distribution function for the third band, fkw(*,3), is not used
 
      data fkw / 0.2747,0.2717,0.2752,0.1177,0.0352,0.0255,
     2           0.1521,0.3974,0.1778,0.1826,0.0374,0.0527,
     3           6*1.00,
     4           0.4654,0.2991,0.1343,0.0646,0.0226,0.0140,
     5           0.5543,0.2723,0.1131,0.0443,0.0160,0.0000,
     6           0.1846,0.2732,0.2353,0.1613,0.1146,0.0310,
     7           0.0740,0.1636,0.4174,0.1783,0.1101,0.0566,
     8           0.1437,0.2197,0.3185,0.2351,0.0647,0.0183/

c-----gkw is the planck-weighted k-distribution function due to h2o
c     line absorption in the 3 subbands (800-720,620-720,540-620 /cm)
c     of band 3 given in table 7.  Note that the order of the sub-bands
c     is reversed.

      data gkw/  0.1782,0.0593,0.0215,0.0068,0.0022,0.0000,
     2           0.0923,0.1675,0.0923,0.0187,0.0178,0.0000,
     3           0.0000,0.1083,0.1581,0.0455,0.0274,0.0041/

 
 
c-----ne is the number of terms used in each band to compute water vapor
c     continuum transmittance (table 6).
 
      data ne /0,0,3,1,1,1,0,0/


c-----tco2 are the six exp factors between levels k1 and k2 
c     tran is the updated total transmittance between levels k1 and k2


c-----th2o is the 6 exp factors between levels k1 and k2 due to
c     h2o line absorption. 

c-----tcon is the 3 exp factors between levels k1 and k2 due to
c     h2o continuum absorption.

c-----trnth2o is the total transmittance between levels k1 and k2 due
c     to both line and continuum absorption computed from eq. (52).


      integer i,j,mmm

      save

ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
         do j=1,n
          do i=1,mmm
           th2o(i,j,1) = th2o(i,j,1)*h2oexp(i,j,k,1)
           th2o(i,j,2) = th2o(i,j,2)*h2oexp(i,j,k,2)
           th2o(i,j,3) = th2o(i,j,3)*h2oexp(i,j,k,3)
           th2o(i,j,4) = th2o(i,j,4)*h2oexp(i,j,k,4)
           th2o(i,j,5) = th2o(i,j,5)*h2oexp(i,j,k,5)
           th2o(i,j,6) = th2o(i,j,6)*h2oexp(i,j,k,6)
          enddo
         enddo


      if (ne(ib).eq.0) then


         do j=1,n
          do i=1,mmm

           trnth2o      =(fkw(1,ib)*th2o(i,j,1)
     *                  + fkw(2,ib)*th2o(i,j,2)
     *                  + fkw(3,ib)*th2o(i,j,3)
     *                  + fkw(4,ib)*th2o(i,j,4)
     *                  + fkw(5,ib)*th2o(i,j,5)
     *                  + fkw(6,ib)*th2o(i,j,6))

          tran(i,j)=tran(i,j)*trnth2o

          enddo
         enddo

      elseif (ne(ib).eq.1) then


         do j=1,n
          do i=1,mmm

           tcon(i,j,1)= tcon(i,j,1)*conexp(i,j,k,1)

           trnth2o      =(fkw(1,ib)*th2o(i,j,1)
     *                  + fkw(2,ib)*th2o(i,j,2)
     *                  + fkw(3,ib)*th2o(i,j,3)
     *                  + fkw(4,ib)*th2o(i,j,4)
     *                  + fkw(5,ib)*th2o(i,j,5)
     *                  + fkw(6,ib)*th2o(i,j,6))*tcon(i,j,1)

          tran(i,j)=tran(i,j)*trnth2o

          enddo
         enddo

      else

         do j=1,n
          do i=1,mmm

           tcon(i,j,1)= tcon(i,j,1)*conexp(i,j,k,1)
           tcon(i,j,2)= tcon(i,j,2)*conexp(i,j,k,2)
           tcon(i,j,3)= tcon(i,j,3)*conexp(i,j,k,3)

           trnth2o      = (  gkw(1,1)*th2o(i,j,1)
     *                     + gkw(2,1)*th2o(i,j,2)
     *                     + gkw(3,1)*th2o(i,j,3)
     *                     + gkw(4,1)*th2o(i,j,4)
     *                     + gkw(5,1)*th2o(i,j,5)
     *                     + gkw(6,1)*th2o(i,j,6) ) * tcon(i,j,1)
     *                  + (  gkw(1,2)*th2o(i,j,1)
     *                     + gkw(2,2)*th2o(i,j,2)
     *                     + gkw(3,2)*th2o(i,j,3)
     *                     + gkw(4,2)*th2o(i,j,4)
     *                     + gkw(5,2)*th2o(i,j,5)
     *                     + gkw(6,2)*th2o(i,j,6) ) * tcon(i,j,2)
     *                  + (  gkw(1,3)*th2o(i,j,1)
     *                     + gkw(2,3)*th2o(i,j,2)
     *                     + gkw(3,3)*th2o(i,j,3)
     *                     + gkw(4,3)*th2o(i,j,4)
     *                     + gkw(5,3)*th2o(i,j,5)
     *                     + gkw(6,3)*th2o(i,j,6) ) * tcon(i,j,3)

          tran(i,j)=tran(i,j)*trnth2o

          enddo
         enddo

      endif



      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine co2kdis(m,n,np,k,co2exp,tco2,tran)
c**********************************************************************
c   compute co2 transmittances between levels k1 and k2 for m x n soundings
c   using the k-distribution method with linear pressure scaling.
c
c   computations follow eq. (34).
c
c---- input parameters
c   number of grid intervals in zonal direction (m)
c   number of grid intervals in meridional direction (n)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2 absorption
c     for the various values of the absorption coefficient (tco2)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,n,np,k

c---- input parameters -----

      real    co2exp(m,n,np,6,2)

c---- updated parameters -----

      real    tco2(m,n,6,2),tran(m,n)

      integer iradave
      COMMON/IPTIONR/ IRADAVE
c---- static data -----

      real    gkc(6,2)

c---- temporary arrays -----

      real    xc

c-----gkc is the planck-weighted co2 k-distribution function 
c     in the band-wing and band-center regions given in table 7.
c     for computing efficiency, sub-bands 3a and 3c are combined.

      data gkc/  0.1395,0.1407,0.1549,0.1357,0.0182,0.0220,
     2           0.0766,0.1372,0.1189,0.0335,0.0169,0.0059/

c-----tco2 is the 6 exp factors between levels k1 and k2. 
c     xc is the total co2 transmittance given by eq. (53).

      integer i,j,mmm

      save

ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
         do j=1,n
          do i=1,mmm

c-----band-wings

           tco2(i,j,1,1)=tco2(i,j,1,1)*co2exp(i,j,k,1,1)
           xc=             gkc(1,1)*tco2(i,j,1,1)

           tco2(i,j,2,1)=tco2(i,j,2,1)*co2exp(i,j,k,2,1)
           xc=xc+gkc(2,1)*tco2(i,j,2,1)

           tco2(i,j,3,1)=tco2(i,j,3,1)*co2exp(i,j,k,3,1)
           xc=xc+gkc(3,1)*tco2(i,j,3,1)

           tco2(i,j,4,1)=tco2(i,j,4,1)*co2exp(i,j,k,4,1)
           xc=xc+gkc(4,1)*tco2(i,j,4,1)

           tco2(i,j,5,1)=tco2(i,j,5,1)*co2exp(i,j,k,5,1)
           xc=xc+gkc(5,1)*tco2(i,j,5,1)

           tco2(i,j,6,1)=tco2(i,j,6,1)*co2exp(i,j,k,6,1)
           xc=xc+gkc(6,1)*tco2(i,j,6,1)

c-----band-center region

           tco2(i,j,1,2)=tco2(i,j,1,2)*co2exp(i,j,k,1,2)
           xc=xc+gkc(1,2)*tco2(i,j,1,2)

           tco2(i,j,2,2)=tco2(i,j,2,2)*co2exp(i,j,k,2,2)
           xc=xc+gkc(2,2)*tco2(i,j,2,2)

           tco2(i,j,3,2)=tco2(i,j,3,2)*co2exp(i,j,k,3,2)
           xc=xc+gkc(3,2)*tco2(i,j,3,2)

           tco2(i,j,4,2)=tco2(i,j,4,2)*co2exp(i,j,k,4,2)
           xc=xc+gkc(4,2)*tco2(i,j,4,2)

           tco2(i,j,5,2)=tco2(i,j,5,2)*co2exp(i,j,k,5,2)
           xc=xc+gkc(5,2)*tco2(i,j,5,2)

           tco2(i,j,6,2)=tco2(i,j,6,2)*co2exp(i,j,k,6,2)
           xc=xc+gkc(6,2)*tco2(i,j,6,2)

           tran(i,j)=tran(i,j)*xc

          enddo
         enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine h2oexps(ib,m,n,np,dh2o,pa,dt,h2oexp)
c**********************************************************************
c   compute exponentials for water vapor line absorption
c   in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer water vapor amount for line absorption (dh2o) 
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  6 exponentials for each layer  (h2oexp)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np

c---- input parameters ------

      real    dh2o(m,n,np),pa(m,n,np),dt(m,n,np)

c---- output parameters -----

      real    h2oexp(m,n,np,6)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- static data -----

      integer mw(8)
      real    xkw(8),aw(8),bw(8)

c---- temporary arrays -----

      real    xh


c-----xkw  are the absorption coefficients for the first
c     k-distribution function due to water vapor line absorption
c     (tables 4 and 7).  units are cm**2/g    
 
      data xkw / 29.55  , 4.167e-1, 1.328e-2, 5.250e-4,
     *            5.25e-4, 2.340e-3, 1.320e-0, 5.250e-4/
 
c-----mw are the ratios between neighboring absorption coefficients
c     for water vapor line absorption (tables 4 and 7).
 
      data mw /6,6,8,6,6,8,6,16/

c-----aw and bw (table 3) are the coefficients for temperature scaling
c     in eq. (25).
 
      data aw/ 0.0021, 0.0140, 0.0167, 0.0302,
     *         0.0307, 0.0154, 0.0008, 0.0096/
      data bw/ -1.01e-5, 5.57e-5, 8.54e-5, 2.96e-4,
     *          2.86e-4, 7.53e-5,-3.52e-6, 1.64e-5/

      integer i,j,k,ik,mmm

      save

c**********************************************************************
c    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
c    and bw.  therefore, h2oexp for these sub-bands are identical.
c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
 
        do k=1,np
         do j=1,n
          do i=1,mmm

c-----xh is   the scaled water vapor amount for line absorption
c     computed from (27).
 
           xh = dh2o(i,j,k)*(pa(i,j,k)*0.002)
     1        * ( 1.+(aw(ib)+bw(ib)* dt(i,j,k))*dt(i,j,k) )

c-----h2oexp is the water vapor transmittance of the layer (k2-1)
c     due to line absorption

           h2oexp(i,j,k,1) = exp(-xh*xkw(ib))

          enddo
         enddo
        enddo

        do ik=2,6

         if(mw(ib).eq.6) then

          do k=1,np
           do j=1,n
            do i=1,mmm
             xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
             h2oexp(i,j,k,ik) = xh*xh*xh
            enddo
           enddo
          enddo

        elseif(mw(ib).eq.8) then

          do k=1,np
           do j=1,n
            do i=1,mmm
             xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
             xh = xh*xh
             h2oexp(i,j,k,ik) = xh*xh
            enddo
           enddo
          enddo

        else

          do k=1,np
           do j=1,n
            do i=1,mmm
             xh = h2oexp(i,j,k,ik-1)*h2oexp(i,j,k,ik-1)
             xh = xh*xh
             xh = xh*xh
             h2oexp(i,j,k,ik) = xh*xh
            enddo
           enddo
          enddo

        endif
       enddo


      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine conexps(ib,m,n,np,dcont,conexp)
c**********************************************************************
c   compute exponentials for continuum absorption in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer scaled water vapor amount for continuum absorption (dcont) 
c
c---- output parameters
c  1 or 3 exponentials for each layer (conexp)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np

c---- input parameters ------

      real    dcont(m,n,np)

c---- updated parameters -----

      real    conexp(m,n,np,3)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- static data -----

      real    xke(8)

c-----xke are the absorption coefficients for the first
c     k-distribution function due to water vapor continuum absorption
c     (table 6).  units are cm**2/g
 
      data xke /  0.00,   0.00,   27.40,   15.8,
     *            9.40,   7.75,     0.0,    0.0/
 
      integer i,j,k,iq,mmm

      save

c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1

        do k=1,np
         do j=1,n
          do i=1,mmm
           conexp(i,j,k,1) = exp(-dcont(i,j,k)*xke(ib))
          enddo
         enddo
        enddo

       if (ib .eq. 3) then

c-----the absorption coefficients for sub-bands 3b (iq=2) and 3a (iq=3)
c     are, respectively, double and quadruple that for sub-band 3c (iq=1)
c     (table 6).

        do iq=2,3
         do k=1,np
          do j=1,n
           do i=1,mmm
            conexp(i,j,k,iq) = conexp(i,j,k,iq-1) *conexp(i,j,k,iq-1)
           enddo
          enddo
         enddo
        enddo

       endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine co2exps(m,n,np,dco2,pa,dt,co2exp)
c       
c**********************************************************************
c   compute co2 exponentials for individual layers.
c
c---- input parameters
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer co2 amount (dco2)
c  layer pressure (pa)
c  layer temperature minus 250K (dt)
c
c---- output parameters
c  6 exponentials for each layer (co2exp)
c**********************************************************************
      implicit none
      integer m,n,np

c---- input parameters -----

      real    dco2(m,n,np),pa(m,n,np),dt(m,n,np)

c---- output parameters -----

      real    co2exp(m,n,np,6,2)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c---- static data -----

      real    xkc(2),ac(2),bc(2),pm(2),prc(2)

c---- temporary arrays -----

      real    xc

c-----xkc is the absorption coefficients for the
c     first k-distribution function due to co2 (table 7).
c     units are 1/(cm-atm)stp.
 
      data xkc/2.656e-5,2.656e-3/
 
c-----parameters (table 3) for computing the scaled co2 amount
c     using (27).

      data prc/  300.0,   30.0/
      data pm /    0.5,   0.85/
      data ac / 0.0182, 0.0042/
      data bc /1.07e-4,2.00e-5/
 
      integer i,j,k,mmm

      save
c**********************************************************************
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
        do k=1,np
         do j=1,n
          do i=1,mmm

c-----compute the scaled co2 amount from eq. (27) for band-wings
c     (sub-bands 3a and 3c).

           xc = dco2(i,j,k)*(pa(i,j,k)/prc(1))**pm(1)
     1             *(1.+(ac(1)+bc(1)*dt(i,j,k))*dt(i,j,k))

c-----six exponential by powers of 8 (table 7).

           co2exp(i,j,k,1,1)=exp(-xc*xkc(1))

           xc=co2exp(i,j,k,1,1)*co2exp(i,j,k,1,1)
           xc=xc*xc
           co2exp(i,j,k,2,1)=xc*xc

           xc=co2exp(i,j,k,2,1)*co2exp(i,j,k,2,1)
           xc=xc*xc
           co2exp(i,j,k,3,1)=xc*xc

           xc=co2exp(i,j,k,3,1)*co2exp(i,j,k,3,1)
           xc=xc*xc
           co2exp(i,j,k,4,1)=xc*xc

           xc=co2exp(i,j,k,4,1)*co2exp(i,j,k,4,1)
           xc=xc*xc
           co2exp(i,j,k,5,1)=xc*xc

           xc=co2exp(i,j,k,5,1)*co2exp(i,j,k,5,1)
           xc=xc*xc
           co2exp(i,j,k,6,1)=xc*xc

c-----compute the scaled co2 amount from eq. (27) for band-center
c     region (sub-band 3b).

           xc = dco2(i,j,k)*(pa(i,j,k)/prc(2))**pm(2)
     1             *(1.+(ac(2)+bc(2)*dt(i,j,k))*dt(i,j,k))

           co2exp(i,j,k,1,2)=exp(-xc*xkc(2))

           xc=co2exp(i,j,k,1,2)*co2exp(i,j,k,1,2)
           xc=xc*xc
           co2exp(i,j,k,2,2)=xc*xc

           xc=co2exp(i,j,k,2,2)*co2exp(i,j,k,2,2)
           xc=xc*xc
           co2exp(i,j,k,3,2)=xc*xc

           xc=co2exp(i,j,k,3,2)*co2exp(i,j,k,3,2)
           xc=xc*xc
           co2exp(i,j,k,4,2)=xc*xc

           xc=co2exp(i,j,k,4,2)*co2exp(i,j,k,4,2)
           xc=xc*xc
           co2exp(i,j,k,5,2)=xc*xc

           xc=co2exp(i,j,k,5,2)*co2exp(i,j,k,5,2)
           xc=xc*xc
           co2exp(i,j,k,6,2)=xc*xc

          enddo
         enddo
        enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c***************  charney  /silo/z1mhy/peng/solar.f ***10/20/95******
      subroutine sorad (m,n,ndim,pl,ta,wa,oa,co2,taucld,reff,taual,
     1             rsirbm,rsuvbm,cosz,flx,fdirir,fdifir,fdirpar,fdifpar)

      implicit none

      integer nxx,nz,nadd,np,nxi,nbb,mm,n1
      PARAMETER (NXX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER (NXI=NXX-2,NBB=1,MM=NXI/NBB/2,N1=1)
      real rflux(mm,7)
      common /radflux/rflux
c********************************************************************
c  
c This routine computes solar fluxes due to the absoption by water
c  vapor, ozone, co2, o2, clouds, and aerosols and due to the
c  scattering by clouds, aerosols, and gases.
c
c This is a vectorized code. It computes the fluxes simultaneous for
c  (m x n) soundings, which is a subset of the (m x ndim) soundings.
c  In a global climate model, m and ndim correspond to the numbers of
c  grid boxes in the zonal and meridional directions, respectively.
c
c Ice and liquid cloud particles are allowed to co-exist in any of the
c  np layers. Two sets of cloud parameters are required as inputs, one
c  for ice paticles and the other for liquid particles.  These parameters
c  are optical thickness (taucld) and effective particle size (reff).
c
c If no information is available for reff, a default value of
c  10 micron for liquid water and 75 micron for ice can be used.
c
c Clouds are grouped into high, middle, and low clouds separated by the
c  level indices ict and icb.  For detail, see the subroutine cldscale.
c
c----- Input parameters:                           
c                                                   units      size
c   number of soundings in zonal direction (m)       n/d        1
c   number of soundings in meridional direction (n)  n/d        1
c   maximum number of soundings in                   n/d        1
c           meridional direction (ndim)
c   number of atmospheric layers (np)                n/d        1
c   level pressure (pl)                              mb       m*ndim*(np+1)
c   layer temperature (ta)                           k        m*ndim*np
c   layer specific humidity (wa)                     gm/gm    m*ndim*np
c   layer ozone concentration (oa)                   gm/gm    m*ndim*np
c   co2 mixing ratio by volumn (co2)               parts/part   1
c   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
c                index 1 for ice particles
c                index 2 for liquid drops
c   effective cloud-particle size (reff)           micrometer m*ndim*np*2
c                index 1 for ice particles
c                index 2 for liquid drops
c   aerosol optical thickness (taual)                n/d      m*ndim*np 
c   solar ir surface albedo for beam                fraction   m*ndim
c                radiation (rsirbm)                
c   uv + par surface albedo for beam                     fraction   m*ndim
c                radiation (rsuvbm)                
c   cosine of solar zenith angle (cosz)            n/d        m*ndim
c
c----- Output parameters
c
c   all-sky flux (downward minus upward) (flx)     fraction   m*ndim*(np+1)
c   all-sky direct downward ir (0.7-10 micron)
c                flux at the surface (fdirir)      fraction   m*ndim
c   all-sky diffuse downward ir flux at
c                the surface (fdifir)              fraction   m*ndim
c   all-sky direct downward par (0.4-0.7 micron)
c                flux at the surface (fdirpar)     fraction   m*ndim
c   all-sky diffuse downward par flux at
c                the surface (fdifpar)             fraction   m*ndim
*
c----- Notes:
c
c    (1) The unit of flux is fraction of the incoming solar radiation
c        at the top of the atmosphere.  Therefore, fluxes should
c        be equal to flux multiplied by the extra-terrestrial solar
c        flux and the cosine of solar zenith angle.
c    (2) Clouds and aerosols can be included in any layers by specifying
c        taucld(i,j,k,*) and taual(i,j,k), k=1,np. 
c        For an atmosphere without clouds and aerosols,
c        set taucld(i,j,k,*)=taual(i,j,k)=0.0.
c    (3) Aerosol single scattering albedos and asymmetry
c        factors are specified in the subroutines solir and soluv.
c    (4) pl(i,j,1) is the pressure at the to of the model, and
c        pl(i,j,np+1) is the surface pressure.
c    (5) the pressure levels ict and icb correspond approximately
c        to 400 and 700 mb.
c        
c**************************************************************************

c-----input parameters

      integer m,n,ndim
      real    pl(m,ndim,np+1),ta(m,ndim,np),wa(m,ndim,np),oa(m,ndim,np)
      real    taucld(m,ndim,np,2),reff(m,ndim,np,2)
      real    taual(m,ndim,np),rsirbm(m,ndim)
      real    rsuvbm(m,ndim),cosz(m,ndim),co2

c-----output parameters

      real    flx(m,ndim,np+1)
      real    fdirir(m,ndim),fdifir(m,ndim)
      real    fdirpar(m,ndim),fdifpar(m,ndim)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----temporary array
 
      integer i,j,k,mmm
      real    dp(MM,N1,np),wh(MM,N1,np),oh(MM,N1,np),scal(MM,N1,np)
      real    swh(MM,N1,np+1),so2(MM,N1,np+1),df(MM,N1,np+1)
      real    sdf(MM,N1),sclr(MM,N1),csm(MM,N1),taux,x,expmn

      save
 
c-----------------------------------------------------------------
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n 
       do i= 1, mmm 

         swh(i,j,1)=0. 
         so2(i,j,1)=0. 

c-----csm is the effective secant of the solar zenith angle
c     see equation (12) of Lacis and Hansen (1974, JAS)    
 
         csm(i,j)=35./sqrt(1224.*cosz(i,j)*cosz(i,j)+1.)

       enddo 
      enddo

      do k= 1, np
       do j= 1, n
         do i= 1, mmm

c-----compute layer thickness and pressure-scaling function. 
c     indices for the surface level and surface layer
c     are np+1 and np, respectively.
 
          dp(i,j,k)=pl(i,j,k+1)-pl(i,j,k)
          scal(i,j,k)=dp(i,j,k)*(.5*(pl(i,j,k)+pl(i,j,k+1))/300.)**.8
 
c-----compute scaled water vapor amount, unit is g/cm**2

          wh(i,j,k)=1.02*wa(i,j,k)*scal(i,j,k)*
     *              (1. + 0.00135*(ta(i,j,k)-240.))
          swh(i,j,k+1)=swh(i,j,k)+wh(i,j,k)

c-----compute ozone amount, unit is (cm-atm)stp.
 
          oh(i,j,k)=1.02*oa(i,j,k)*dp(i,j,k)*466.7

        enddo
       enddo
      enddo

c-----initialize fluxes for all-sky (flx) and flux reduction (df)

      do k=1, np+1
       do j=1, n
        do i=1, mmm
          flx(i,j,k)=0.
          df(i,j,k)=0.
        enddo
       enddo
      enddo

c-----compute solar ir fluxes

      call solir (m,n,ndim,wh,taucld,reff,taual,
     *            csm,rsirbm,flx,fdirir,fdifir)

c-----compute and update uv and par fluxes

      call soluv (m,n,ndim,oh,dp,taucld,reff,taual,
     *            csm,rsuvbm,flx,fdirpar,fdifpar)

c-----compute scaled amount of o2 (so2), unit is (cm-atm)stp.

      do k= 1, np
       do j= 1, n
        do i= 1, mmm
          so2(i,j,k+1)=so2(i,j,k)+165.22*scal(i,j,k)
        enddo
       enddo
      enddo

c-----compute flux reduction due to oxygen following
c      chou (J. climate, 1990). The fraction 0.0287 is the
c      extraterrestrial solar flux in the o2 bands.

       do k= 2, np+1
        do j= 1, n
         do i= 1, mmm
           x=so2(i,j,k)*csm(i,j)
           df(i,j,k)=df(i,j,k)+0.0287*(1.-expmn(-0.00027*sqrt(x)))
         enddo
        enddo
       enddo          

c-----compute scaled amounts for co2 (so2). unit is (cm-atm)stp.

      do k= 1, np
       do j= 1, n
        do i= 1, mmm
         so2(i,j,k+1)=so2(i,j,k)+co2*789.*scal(i,j,k)
        enddo
       enddo
      enddo

c-----compute and update flux reduction due to co2 following
c     chou (J. Climate, 1990)

      call flxco2(m,n,np,so2,swh,csm,df)

c-----adjust for the all-sky fluxes due to o2 and co2.  It is
c     assumed that o2 and co2 have no effects on solar radiation
c     below clouds.

      do j=1,n
       do i=1,mmm
          sdf(i,j)=0.0 
          sclr(i,j)=1.0 
       enddo
      enddo

      do k=1,np
       do j=1,n
        do i=1,mmm

           taux=taucld(i,j,k,1)+taucld(i,j,k,2)
         if(taux.gt.0.01 .and. sclr(i,j).eq.1.) then
          sdf(i,j)=df(i,j,k)
          sclr(i,j)=0.0
         endif

          flx(i,j,k+1)=flx(i,j,k+1)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)
 
C --------------- GCSS -------------------------------------------
         if(k.eq.1)then
           rflux(i,3)=rflux(i,3)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)
         endif
         if(k.eq.(np+1))then
           rflux(i,1)=rflux(i,1)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)
         endif
C ---------------------------------------------------------------
         


        enddo
       enddo
      enddo

c-----adjust for the direct downward ir flux.
      do j= 1, n
       do i= 1, mmm
           fdirir(i,j)=fdirir(i,j)-sdf(i,j)-df(i,j,np+1)*sclr(i,j)
       enddo
      enddo

      return
      end  


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solir (m,n,ndim,wh,taucld,reff,taual,
     $                  csm,rsirbm,flx,fdirir,fdifir)

      implicit none

      integer nxx,nz,nadd,np,nxi,nbb,mm,n1,nk,nband
      PARAMETER(NXX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER(NXI=NXX-2)
      PARAMETER(NBB=1,MM=NXI/NBB/2)
      PARAMETER(N1=1)
      parameter (nk=10,nband=3)
      real rflux(mm,7)
      common/radflux/rflux
c************************************************************************
c  compute solar flux in the infrared region. The spectrum is divided
c   into three bands:
c
c          band   wavenumber(/cm)  wavelength (micron)
c           1       1000-4400         2.27-10.0
c           2       4400-8200         1.22-2.27
c           3       8200-14300        0.70-1.22
c
c----- Input parameters:                            units      size
c
c   number of soundings in zonal direction (m)       n/d        1
c   number of soundings in meridional direction (n)  n/d        1
c   maximum number of soundings in                   n/d        1
c          meridional direction (ndim)
c   number of atmospheric layers (np)                n/d        1
c   layer water vapor content (wh)                 gm/cm^2    m*n*np
c   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   effective cloud-particle size (reff)           micrometer m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   aerosol optical thickness (taual)                n/d      m*ndim*np
c   cosecant of the solar zenith angle (csm)         n/d      m*n
c   near ir surface albedo for beam                fraction   m*ndim
c                radiation (rsirbm)
c
c----- output (updated) parameters:
c
c   all-sky flux (downward-upward) (flx)           fraction   m*ndim*(np+1)
c   all-sky direct downward ir flux at
c          the surface (fdirir)                    fraction   m*ndim
c   all-sky diffuse downward ir flux at
c          the surface (fdifir)                    fraction   m*ndim
c
c----- note: the following parameters must be specified by users:
c   aerosol single scattering albedo (ssaal)         n/d      nband
c   aerosol asymmetry factor (asyal)                 n/d      nband
c
c*************************************************************************

c-----input parameters

      integer m,n,ndim
      real    taucld(m,ndim,np,2),reff(m,ndim,np,2),rsirbm(m,ndim)
      real    wh(m,n,np),taual(m,ndim,np),csm(m,n)

c-----output (updated) parameters

      real    flx(m,ndim,np+1)
      real    fdirir(m,ndim),fdifir(m,ndim)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----static parameters

      real    xk(nk),hk(nband,nk),ssaal(nband),asyal(nband)
      real    aia(nband,3),awa(nband,3),aig(nband,3),awg(nband,3)

c-----temporary array

      integer ib,ik,i,j,k,mmm
      real    ssacl(MM,N1,np),asycl(MM,N1,np)
      real    rr(MM,N1,np+1),tt(MM,N1,np+1),td(MM,N1,np+1),
     $        rs(MM,N1,np+1),ts(MM,N1,np+1)
      real    flxdn(MM,N1,np+1),fdndir(MM,N1),fdndif(MM,N1)
      real    tauwv,tausto,ssatau,asysto,tauto,ssato,asyto
      real    taux,reff1,reff2,w1,w2,g1,g2
      real    ssaclt(MM,N1),asyclt(MM,N1)
      real   rr1t(MM,N1),tt1t(MM,N1),td1t(MM,N1),rs1t(MM,N1),ts1t(MM,N1)

c-----water vapor absorption coefficient for 10 k-intervals.
c     unit: cm^2/gm

      data xk/            
     1  0.0010, 0.0133, 0.0422, 0.1334, 0.4217,            
     2  1.334,  5.623,  31.62,  177.8,  1000.0/  

c-----water vapor k-distribution function,
c     the sum of hk is 0.52926. unit: fraction

      data hk/
     1 .01074,.08236,.20673,  .00360,.01157,.03497,
     2 .00411,.01133,.03011,  .00421,.01143,.02260,
     3 .00389,.01240,.01336,  .00326,.01258,.00696,
     4 .00499,.01381,.00441,  .00465,.00650,.00115,
     5 .00245,.00244,.00026,  .00145,.00094,.00000/

c-----aerosol single-scattering albedo and asymmetry factor

      data ssaal/0.999, 0.999, 0.999/
      data asyal/0.850, 0.850, 0.850/
 
c-----coefficients for computing the single scattering albedo of
c     ice clouds from ssa=1-(aia(*,1)+aia(*,2)*reff+aia(*,3)*reff**2)

      data aia/
     1  .08938331, .00215346,-.00000260,
     2  .00299387, .00073709, .00000746,
     3 -.00001038,-.00000134, .00000000/

c-----coefficients for computing the single scattering albedo of
c     liquid clouds from ssa=1-(awa(*,1)+awa(*,2)*reff+awa(*,3)*reff**2)

      data awa/
     1  .01209318,-.00019934, .00000007,
     2  .01784739, .00088757, .00000845,
     3 -.00036910,-.00000650,-.00000004/

c-----coefficients for computing the asymmetry factor of ice clouds
c     from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

      data aig/
     1  .84090400, .76098937, .74935228,
     2  .00126222, .00141864, .00119715,
     3 -.00000385,-.00000396,-.00000367/

c-----coefficients for computing the asymmetry factor of liquid clouds
c     from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

      data awg/
     1  .83530748, .74513197, .79375035,
     2  .00257181, .01370071, .00832441,
     3  .00005519,-.00038203,-.00023263/

      save

c-----initialize surface fluxes, reflectances, and transmittances
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n
       do i= 1, mmm
         fdirir(i,j)=0.0
         fdifir(i,j)=0.0
         rr(i,j,np+1)=rsirbm(i,j)
         rs(i,j,np+1)=rsirbm(i,j)
         td(i,j,np+1)=0.0
         tt(i,j,np+1)=0.0
         ts(i,j,np+1)=0.0
       enddo
      enddo

c-----integration over spectral bands

      do 100 ib=1,nband

c-----compute cloud single scattering albedo and asymmetry factor
c     for a mixture of ice and liquid particles.

       do k= 1, np

        do j= 1, n
         do i= 1, mmm

           ssaclt(i,j)=1.0
           asyclt(i,j)=1.0

           taux=taucld(i,j,k,1)+taucld(i,j,k,2)
          if (taux.gt.0.05) then

           reff1=min(reff(i,j,k,1),130.)
           reff2=min(reff(i,j,k,2),20.0)

           w1=(1.-(aia(ib,1)+(aia(ib,2)+
     *         aia(ib,3)*reff1)*reff1))*taucld(i,j,k,1)
           w2=(1.-(awa(ib,1)+(awa(ib,2)+
     *         awa(ib,3)*reff2)*reff2))*taucld(i,j,k,2)
           ssaclt(i,j)=(w1+w2)/taux

           g1=(aig(ib,1)+(aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
           g2=(awg(ib,1)+(awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
           asyclt(i,j)=(g1+g2)/(w1+w2)

          endif

         enddo
        enddo

        do j=1,n
         do i=1,mmm
            ssacl(i,j,k)=ssaclt(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            asycl(i,j,k)=asyclt(i,j)
         enddo
        enddo

       enddo

c-----integration over the k-distribution function

         do 200 ik=1,nk

          do 300 k= 1, np

           do j= 1, n
            do i= 1, mmm

             tauwv=xk(ik)*wh(i,j,k)
 
c-----compute total optical thickness, single scattering albedo,
c     and asymmetry factor.
 
             tausto=tauwv+taual(i,j,k)+1.0e-8
             ssatau=ssaal(ib)*taual(i,j,k)
             asysto=asyal(ib)*ssaal(ib)*taual(i,j,k)
 
c-----compute reflectance and transmittance

              taux=taucld(i,j,k,1)+taucld(i,j,k,2)
              tauto=tausto+taux
c             tauto=max(tausto+taux,1.0e-8)
              ssato=(ssatau+ssacl(i,j,k)*taux)/tauto+1.0e-8
c             ssato=max( (ssatau+ssacl(i,j,k)*taux)/tauto,1.0e-8)
              ssato=min(ssato,0.999999)
              asyto=(asysto+asycl(i,j,k)*ssacl(i,j,k)*taux)/
     *              (ssato*tauto)

              call deledd (tauto,ssato,asyto,csm(i,j), 
     *                     rr1t(i,j),tt1t(i,j),td1t(i,j))

              call sagpol (tauto,ssato,asyto,rs1t(i,j),ts1t(i,j))

            enddo
           enddo

           do j=1,n
            do i=1,mmm
               rr(i,j,k)=rr1t(i,j)
            enddo
           enddo
           do j=1,n
            do i=1,mmm
               tt(i,j,k)=tt1t(i,j)
            enddo
           enddo
           do j=1,n
            do i=1,mmm
               td(i,j,k)=td1t(i,j)
            enddo
           enddo
           do j=1,n
            do i=1,mmm
               rs(i,j,k)=rs1t(i,j)
            enddo
           enddo
           do j=1,n
            do i=1,mmm
               ts(i,j,k)=ts1t(i,j)
            enddo
           enddo

 300  continue

c-----flux calculations at each level using the two-stream adding method
 
       call adding (m,n,rr,tt,td,rs,ts,flxdn,fdndir,fdndif)

       do k= 1, np+1
        do j= 1, n
         do i= 1, mmm
          flx(i,j,k) = flx(i,j,k)+flxdn(i,j,k)*hk(ib,ik)

C --------- GCSS ---------------------------------------
         if(k.eq.1)then
           rflux(i,3)=rflux(i,3)+hk(ib,ik)
         endif
         if(k.eq.(np+1))then
           rflux(i,1)=rflux(i,1)+rflux(i,2)*hk(ib,ik)
         endif
C -------------------------------------------------------

         enddo
        enddo
       enddo

       do j= 1, n
        do i= 1, mmm
          fdirir(i,j) = fdirir(i,j)+fdndir(i,j)*hk(ib,ik)
          fdifir(i,j) = fdifir(i,j)+fdndif(i,j)*hk(ib,ik)
        enddo
       enddo

  200 continue
  100 continue
 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine soluv (m,n,ndim,oh,dp,taucld,reff,taual,
     *                  csm,rsuvbm,flx,fdirpar,fdifpar)

      implicit none

      integer nxx,nz,nadd,np,nxi,nbb,mm,n1,nband
      PARAMETER(NXX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER(NXI=NXX-2,NBB=1,MM=NXI/NBB/2,N1=1)
      parameter (nband=8)
      real rflux(mm,7)
      common/radflux/rflux
c************************************************************************
c  compute solar fluxes in the uv+visible region. the spectrum is
c  grouped into 8 bands:
c  
c              Band     Micrometer
c
c       UV-C    1.     .175 - .225
c               2.     .225 - .245
c                      .260 - .280
c               3.     .245 - .260
c
c       UV-B    4.     .280 - .295
c               5.     .295 - .310
c               6.     .310 - .320
c      
c       UV-A    7.     .320 - .400
c      
c       PAR     8.     .400 - .700
c
c----- Input parameters:                            units      size
c
c   number of soundings in zonal direction (m)       n/d        1
c   number of soundings in meridional direction (n)  n/d        1
c   maximum number of soundings in                   n/d        1
c           meridional direction (ndim)
c   number of atmospheric layers (np)                n/d        1
c   layer ozone content (oh)                      (cm-atm)stp m*n*np
c   layer pressure thickness (dp)                    mb       m*n*np
c   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   effective cloud-particle size (reff)           micrometer m*ndim*np*2
c          index 1 for ice paticles
c          index 2 for liquid particles
c   aerosol optical thickness (taual)                n/d      m*ndim*np
c   cosecant of the solar zenith angle (csm)         n/d      m*n
c   uv+par surface albedo for beam                 fraction   m*ndim
c           radiation (rsuvbm)
c
c----- output (updated) parameters:
c
c   all-sky net downward flux (flx)                fraction   m*ndim*(np+1)
c   all-sky direct downward par flux at
c          the surface (fdirpar)                   fraction   m*ndim
c   all-sky diffuse downward par flux at
c          the surface (fdifpar)                   fraction   m*ndim
c
c----- note: the following parameters must be specified by users:
c
c   aerosol single scattering albedo (ssaal)         n/d        1
c   aerosol asymmetry factor (asyal)                 n/d        1
c
*
c***********************************************************************

c-----input parameters

      integer m,n,ndim
      real    taucld(m,ndim,np,2),reff(m,ndim,np,2)
      real    oh(m,n,np),dp(m,n,np),taual(m,ndim,np)
      real    rsuvbm(m,ndim),csm(m,n)

c-----output (updated) parameter

      real    flx(m,ndim,np+1)
      real    fdirpar(m,ndim),fdifpar(m,ndim)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----static parameters

c      integer nband

      real    hk(nband),xk(nband),ry(nband)
      real    asyal(nband),ssaal(nband),aig(3),awg(3)

c-----temporary array

      integer i,j,k,ib,mmm
      real    taurs,tauoz,tausto,ssatau,asysto,tauto,ssato,asyto
      real    taux,reff1,reff2,g1,g2,asycl(MM,N1,np)
      real    td(MM,N1,np+1),rr(MM,N1,np+1),tt(MM,N1,np+1),
     $        rs(MM,N1,np+1),ts(MM,N1,np+1)
      real    flxdn(MM,N1,np+1),fdndir(MM,N1),fdndif(MM,N1)
      real    asyclt(MM,N1)
      real   rr1t(MM,N1),tt1t(MM,N1),td1t(MM,N1),rs1t(MM,N1),ts1t(MM,N1)

c-----hk is the fractional extra-terrestrial solar flux.
c     the sum of hk is 0.47074.

      data hk/.00057, .00367, .00083, .00417,
     *        .00600, .00556, .05913, .39081/

c-----xk is the ozone absorption coefficient. unit: /(cm-atm)stp

      data xk /30.47, 187.2,  301.9,   42.83,
     *         7.09,  1.25,   0.0345,  0.0539/

c-----ry is the extinction coefficient for Rayleigh scattering.
c     unit: /mb.

      data ry /.00604, .00170, .00222, .00132,
     *         .00107, .00091, .00055, .00012/

c-----aerosol single-scattering albedo and asymmetry factor

      data ssaal/0.999,0.999,0.999,0.999,0.999,0.999,0.999,0.999/
      data asyal/0.850,0.850,0.850,0.850,0.850,0.850,0.850,0.850/

c-----coefficients for computing the asymmetry factor of ice clouds
c     from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

      data aig/.74625000,.00105410,-.00000264/

c-----coefficients for computing the asymmetry factor of liquid
c     clouds from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

      data awg/.82562000,.00529000,-.00014866/

      save

c-----initialize surface reflectances and transmittances
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n
       do i= 1, mmm                    
         rr(i,j,np+1)=rsuvbm(i,j)
         rs(i,j,np+1)=rsuvbm(i,j)
         td(i,j,np+1)=0.0
         tt(i,j,np+1)=0.0
         ts(i,j,np+1)=0.0
       enddo
      enddo

c-----compute cloud asymmetry factor for a mixture of
c     liquid and ice particles.  unit of reff is micrometers.

      do k= 1, np

       do j= 1, n
        do i= 1, mmm

           asyclt(i,j)=1.0

           taux=taucld(i,j,k,1)+taucld(i,j,k,2)
          if (taux.gt.0.05) then

           reff1=min(reff(i,j,k,1),130.)
           reff2=min(reff(i,j,k,2),20.0)

           g1=(aig(1)+(aig(2)+aig(3)*reff1)*reff1)*taucld(i,j,k,1)
           g2=(awg(1)+(awg(2)+awg(3)*reff2)*reff2)*taucld(i,j,k,2)
           asyclt(i,j)=(g1+g2)/taux

          endif

        enddo
       enddo

       do j=1,n
        do i=1,mmm
           asycl(i,j,k)=asyclt(i,j)
        enddo
       enddo

      enddo
            
c-----integration over spectral bands

      do 100 ib=1,nband

       do 300 k= 1, np

        do j= 1, n
         do i= 1, mmm

c-----compute ozone and rayleigh optical thicknesses

          taurs=ry(ib)*dp(i,j,k)
          tauoz=xk(ib)*oh(i,j,k)
 
c-----compute total optical thickness, single scattering albedo,
c     and asymmetry factor

          tausto=taurs+tauoz+taual(i,j,k)+1.0e-8
          ssatau=ssaal(ib)*taual(i,j,k)+taurs
          asysto=asyal(ib)*ssaal(ib)*taual(i,j,k)

c-----compute reflectance and transmittance

            taux=taucld(i,j,k,1)+taucld(i,j,k,2)
           tauto=tausto+taux
           ssato=(ssatau+taux)/tauto+1.0e-8
           ssato=min(ssato,0.999999)
           asyto=(asysto+asycl(i,j,k)*taux)/(ssato*tauto)

           call deledd (tauto,ssato,asyto,csm(i,j), 
     *                  rr1t(i,j),tt1t(i,j),td1t(i,j))

           call sagpol (tauto,ssato,asyto,rs1t(i,j),ts1t(i,j))

         enddo
        enddo

        do j=1,n
         do i=1,mmm
            rr(i,j,k)=rr1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            tt(i,j,k)=tt1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            td(i,j,k)=td1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            rs(i,j,k)=rs1t(i,j)
         enddo
        enddo
        do j=1,n
         do i=1,mmm
            ts(i,j,k)=ts1t(i,j)
         enddo
        enddo

 300  continue

c-----flux calculations

       call adding (m,n,rr,tt,td,rs,ts,flxdn,fdndir,fdndif)

       do k= 1, np+1
        do j= 1, n
         do i= 1, mmm
          flx(i,j,k)=flx(i,j,k)+flxdn(i,j,k)*hk(ib)

C --------- GCSS ---------------------------------------
         if(k.eq.1)then
           rflux(i,3)=rflux(i,3)+hk(ib)
         endif
         if(k.eq.(np+1))then
           rflux(i,1)=rflux(i,1)+rflux(i,2)*hk(ib)
         endif
C -------------------------------------------------------

         enddo
        enddo
       enddo

       if(ib.eq.8) then
         do j=1,n
          do i=1,mmm
           fdirpar(i,j)=fdndir(i,j)*hk(ib)
           fdifpar(i,j)=fdndif(i,j)*hk(ib)
         enddo
        enddo
       endif

 100  continue

      return
      end

c*********************************************************************

      subroutine deledd(tau,ssc,g0,csm,rr,tt,td)

c*********************************************************************
c
c-----uses the delta-eddington approximation to compute the
c     bulk scattering properties of a single layer
c     coded following King and Harshvardhan (JAS, 1986)
c
c  inputs:
c
c     tau: the effective optical thickness
c     ssc: the effective single scattering albedo
c     g0:  the effective asymmetry factor
c     csm: the effective secant of the zenith angle
c
c  outputs:
c
c     rr: the layer reflection of the direct beam
c     tt: the layer diffuse transmission of the direct beam
c     td: the layer direct transmission of the direct beam
c
c*********************************************************************

      implicit none

      real zero,one,two,three,four,fourth,seven,thresh
      parameter (one =1., three=3.)
      parameter (two =2., seven=7.)
      parameter (four=4., fourth=.25)
      parameter (zero=0., thresh=1.e-8)
      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----input parameters
      real tau,ssc,g0,csm

c-----output parameters
      real rr,tt,td

c-----temporary parameters

      real zth,ff,xx,taup,sscp,gp,gm1,gm2,gm3,akk,alf1,alf2,
     *     all,bll,st7,st8,cll,dll,fll,ell,st1,st2,st3,st4
      save

c---------------------------------------------------------------------

                zth = one / csm

c  delta-eddington scaling of single scattering albedo,
c  optical thickness, and asymmetry factor,
c  K & H eqs(27-29)

                ff  = g0*g0
                xx  = one-ff*ssc
                taup= tau*xx
                sscp= ssc*(one-ff)/xx
                gp  = g0/(one+g0)

c  gamma1, gamma2, and gamma3. see table 2 and eq(26) K & H
c  ssc and gp are the d-s single scattering
c  albedo and asymmetry factor.

                xx  =  three*gp
                gm1 =  (seven - sscp*(four+xx))*fourth
                gm2 = -(one   - sscp*(four-xx))*fourth

c  akk is k as defined in eq(25) of K & H

                akk = sqrt((gm1+gm2)*(gm1-gm2))

                xx  = akk * zth
                st7 = one - xx
                st8 = one + xx
                st3 = st7 * st8

                if (abs(st3) .lt. thresh) then
                    zth = zth + 0.001
                    xx  = akk * zth
                    st7 = one - xx
                    st8 = one + xx
                    st3 = st7 * st8
                endif

c  extinction of the direct beam transmission

                td  = exp(-taup/zth)

c  alf1 and alf2 are alpha1 and alpha2 from eqs (23) & (24) of K & H

                gm3  = (two - zth*three*gp)*fourth
                xx   = gm1 - gm2
                alf1 = gm1 - gm3 * xx
                alf2 = gm2 + gm3 * xx

c  all is last term in eq(21) of K & H
c  bll is last term in eq(22) of K & H

                xx  = akk * two
                all = (gm3 - alf2 * zth    )*xx*td
                bll = (one - gm3 + alf1*zth)*xx

                xx  = akk * gm3
                cll = (alf2 + xx) * st7
                dll = (alf2 - xx) * st8

                xx  = akk * (one-gm3)
                fll = (alf1 + xx) * st8
                ell = (alf1 - xx) * st7

                st2 = exp(-akk*taup)
                st4 = st2 * st2

                st1 =  sscp / ((akk+gm1 + (akk-gm1)*st4) * st3)

c  rr is r-hat of eq(21) of K & H
c  tt is diffuse part of t-hat of eq(22) of K & H

                rr =   ( cll-dll*st4    -all*st2)*st1
                tt = - ((fll-ell*st4)*td-bll*st2)*st1

                rr = max(rr,zero)
                tt = max(tt,zero)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sagpol(tau,ssc,g0,rll,tll)

c*********************************************************************
c-----transmittance (tll) and reflectance (rll) of diffuse radiation
c     follows Sagan and Pollock (JGR, 1967).
c     also, eq.(31) of Lacis and Hansen (JAS, 1974).
c
c-----input parameters:
c
c      tau: the effective optical thickness
c      ssc: the effective single scattering albedo
c      g0:  the effective asymmetry factor
c
c-----output parameters:
c
c      rll: the layer reflection of diffuse radiation
c      tll: the layer transmission of diffuse radiation
c
c*********************************************************************

      implicit none

      real    one,three,four
      parameter (one=1., three=3., four=4.)

c-----output parameters:

      real    tau,ssc,g0,expmn

c-----output parameters:

      real    rll,tll

c-----temporary arrays

      real    xx,uuu,ttt,emt,up1,um1,st1

      save

          if(ssc .gt. 0.001) then

             xx  = one-ssc*g0
             uuu = sqrt( xx/(one-ssc))
             ttt = sqrt( xx*(one-ssc)*three )*tau
             emt = expmn(-ttt)
             up1 = uuu + one
             um1 = uuu - one
             xx  = um1*emt
             st1 = one / ((up1+xx) * (up1-xx))
             rll = up1*um1*(one-emt*emt)*st1
             tll = uuu*four*emt         *st1

          else

             rll = 0.0
             tll = expmn(-1.66*tau)

          endif

      return
      end

c*******************************************************************

      function expmn(fin)

      implicit none
c*******************************************************************
c compute exponential for arguments in the range 0> fin > -10.

      real    one,expmin,e1,e2,e3,e4
      parameter (one=1.0, expmin=-10.0)
      parameter (e1=1.0,        e2=-2.507213e-1)
      parameter (e3=2.92732e-2, e4=-3.827800e-3)

      real    fin,expmn

      if (fin .lt. expmin) fin = expmin
      expmn = ((e4*fin + e3)*fin+e2)*fin+e1
      expmn = expmn * expmn
      expmn = one / (expmn * expmn)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine adding (m,n,rr,tt,td,rs,ts,flxdn,fdndir,fdndif)

      implicit none

      integer nx,nz,nadd,np,nxi,nbb,mm,n1
      PARAMETER(NX=514,NZ=43,NADD=7,NP=NZ-2+NADD)
      PARAMETER(NXI=NX-2)
      PARAMETER(NBB=1,MM=NXI/NBB/2)
      PARAMETER(N1=1)
      real rflux(mm,7)
      common/radflux/rflux
c*********************************************************************
c  compute upward and downward fluxes using a two-stream adding method
c  computations follow equations (3)-(5) of Chou (1992, JAS)
c  input parameters:
c     np: total number of layers 
c     rr:  reflection of a layer illuminated by beam radiation
c     tt:  diffuse transmission of a layer illuminated by beam radiation
c     td: direct beam tranmssion
c     ts: transmission of a layer illuminated by diffuse radiation
c     rs: reflection of a layer illuminated by diffuse radiation
c
c  output parameters:
c     flxdn:  net downward fluxes
c     fdndir: surface direct downward flux
c     fdndif: surface diffuse downward flux
c*********************************************************************

c-----input parameters

      integer m,n
      real    rr(m,n,np+1),tt(m,n,np+1),td(m,n,np+1),rs(m,n,np+1),
     $        ts(m,n,np+1)

c-----output parameters

      real    flxdn(m,n,np+1),fdndir(m,n),fdndif(m,n)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----temporary array

      integer i,j,k,mmm
      real    denm,xx,fupdif
      real    rssab(MM,N1,np+1),rabx(MM,N1,np+1),rsabx(MM,N1,np+1)
      real    tbab(MM,N1,np+1),tab(MM,N1,np+1)

      save

c-----layers are added one at a time, going down from the top layer,
c     tbab is the composite transmittance illuminated by beam radiation
c     tab is the composite diffuse transmittance illuminated by
c         beam radiation
c     rssab is the composite reflectance illuminated from below
c         by diffuse radiation
c     tab and rssab are computed from eqs. (4b) and (3b) of Chou
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do j= 1, n
        do i= 1, mmm
          tbab(i,j,1)  = td(i,j,1)
          tab(i,j,1)   = tt(i,j,1)
          rssab(i,j,1) = rs(i,j,1)
        enddo
      enddo

      do k= 2, np
       do j= 1, n
        do i= 1, mmm
                 denm = ts(i,j,k)/( 1.-rssab(i,j,k-1)*rs(i,j,k) )
          tbab(i,j,k) = tbab(i,j,k-1)*td(i,j,k)
          tab(i,j,k)  = tbab(i,j,k-1)*tt(i,j,k)
     1                + (tbab(i,j,k-1)*rssab(i,j,k-1)
     2                * rr(i,j,k)+tab(i,j,k-1))*denm
          rssab(i,j,k)= rs(i,j,k)+ts(i,j,k)*rssab(i,j,k-1)*denm
        enddo
       enddo
      enddo

c-----layers are added one at a time, going up
c     rabx is the composite reflectance illuminated by beam radiation
c     rsabx is the composite reflectance illuminated from above
c         by diffuse radiation
c     rabx and rsabx are computed from eqs. (4a) and (3a) of Chou
 
      do j= 1, n
       do i= 1, mmm
         rabx(i,j,np+1) = rr(i,j,np+1)
         rsabx(i,j,np+1)= rs(i,j,np+1)
       enddo
      enddo

      do k= np, 1, -1
       do j= 1, n
        do i= 1, mmm
                 denm  = ts(i,j,k)/( 1.-rs(i,j,k)*rsabx(i,j,k+1) )
          rabx(i,j,k)  = rr(i,j,k)+(td(i,j,k)*rabx(i,j,k+1)
     *                 + tt(i,j,k)*rsabx(i,j,k+1))*denm
          rsabx(i,j,k) = rs(i,j,k)+ts(i,j,k)*rsabx(i,j,k+1)*denm
        enddo
       enddo
      enddo
 
c-----compute fluxes following eq (5) of Chou (1992)
 
c     fdndir is the direct  downward flux
c     fdndif is the diffuse downward flux
c     fupdif is the diffuse upward flux

      do k=2,np+1
       do j=1, n
        do i=1, mmm
                denm  = 1./(1.- rssab(i,j,k-1)*rsabx(i,j,k))
         fdndir(i,j)  = tbab(i,j,k-1)
                  xx  = tbab(i,j,k-1)*rabx(i,j,k)
         fdndif(i,j)  = (xx*rssab(i,j,k-1)+tab(i,j,k-1))*denm
              fupdif  = (xx+tab(i,j,k-1)*rsabx(i,j,k))*denm
         flxdn(i,j,k) = fdndir(i,j)+fdndif(i,j)-fupdif

C ------------- GCSS ---------------------------------------------------
         if(k.eq.(np+1))then
           rflux(i,2)=fdndir(i,j)+fdndif(i,j)
         endif
C ----------------------------------------------------------------------
         
        enddo
       enddo
      enddo

       do j=1, n
        do i=1, mmm
         flxdn(i,j,1) = 1.0-rabx(i,j,1)
        enddo
       enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine flxco2(m,n,np,swc,swh,csm,df)

c*****************************************************************

c-----compute the reduction of clear-sky downward solar flux
c     due to co2 absorption.

      implicit none

c-----input parameters

      integer m,n,np
      real    csm(m,n),swc(m,n,np+1),swh(m,n,np+1),cah(22,19)

c-----output (undated) parameter

      real    df(m,n,np+1)

      integer iradave
      COMMON/IPTIONR/ IRADAVE

c-----temporary array

      integer i,j,k,ic,iw,mmm
      real    xx,clog,wlog,dc,dw,x1,x2,y2

c********************************************************************
c-----include co2 look-up table

      include "cah.dat"

      save
 
c********************************************************************
c-----table look-up for the reduction of clear-sky solar
c     radiation due to co2. The fraction 0.0343 is the
c     extraterrestrial solar flux in the co2 bands.
ctao
      MMM=M
      IF (IRADAVE .EQ. 1) MMM=1
ctao
      do k= 2, np+1
       do j= 1, n
        do i= 1, mmm
          xx=1./.3
          clog=log10(swc(i,j,k)*csm(i,j))
          wlog=log10(swh(i,j,k)*csm(i,j))
          ic=int( (clog+3.15)*xx+1.)
          iw=int( (wlog+4.15)*xx+1.)
          if(ic.lt.2)ic=2
          if(iw.lt.2)iw=2
          if(ic.gt.22)ic=22
          if(iw.gt.19)iw=19     
          dc=clog-float(ic-2)*.3+3.
          dw=wlog-float(iw-2)*.3+4.   
          x1=cah(1,iw-1)+(cah(1,iw)-cah(1,iw-1))*xx*dw
          x2=cah(ic-1,iw-1)+(cah(ic-1,iw)-cah(ic-1,iw-1))*xx*dw
          y2=x2+(cah(ic,iw-1)-cah(ic-1,iw-1))*xx*dc
          if (x1.lt.y2) x1=y2
          df(i,j,k)=df(i,j,k)+0.0343*(x1-y2)
        enddo     
       enddo      
      enddo      

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine budget (taa)

      implicit none
      integer nx,nz,itt
      parameter (NX=514,NZ=43,ITT=244)

      real    taa(nz)

      integer imax,iles,il2,kmax,kles,kl2
      common/bxz/ imax,iles,il2,kmax,kles,kl2

      real    dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2

      real    dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2

      real    tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b5/ tb,qb,rho1,rho,ba,bb,ta,
     $           qa,ta1,qa1,coef,c1,c2,c3,am,
     $           am1,ub,vb,wb,ub1,vb1,rrho,rrho1
      real    tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/b6/ tls,qls,fd,fe,p0,pi,f0,st,sv,sq,sc,se,sqa,ri,ar,rx

      real    smf(nz),smu(nz),smd(nz),smn(nz),stf(nz),stu(nz),
     $   std(nz),stn(nz),sqf(nz),squ(nz),sqd(nz),sqn(nz),sdt(nz),
     $   sdq(nz),stv(nz),
     $   stf1(nz),stu1(nz),std1(nz),stn1(nz),sqf1(nz),squ1(nz),
     $   sqd1(nz),sqn1(nz),stf2(nz),stu2(nz),std2(nz),stn2(nz),sqf2(nz),
     $   squ2(nz),sqd2(nz),sqn2(nz),cld(nz),cldu(nz),cldd(nz),cldn(nz)
      common/b8/ smf,smu,smd,smn,stf,stu,std,stn,sqf,squ,sqd,sqn,sdt,
     $           sdq,stv,stf1,stu1,std1,stn1,sqf1,squ1,
     $           sqd1,sqn1,stf2,stu2,std2,stn2,sqf2,
     $           squ2,sqd2,sqn2,cld,cldu,cldd,cldn

      real    tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     1  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     2  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     3  q2t(nz),vb_6h(nz,itt),vbt(nz)
      common/bb6/ tls1,tls2,qls1,qls2,tls3,tls4,qls3,qls4,sft,sfq,wbt,
     $            wb_6h,ub_6h,ubt,q1_6h,q1t,q2_6h,q2t,vb_6h,vbt
      real    thom(nz,4,7),tdw(nz,4,7),tmlt(nz,4,7),saut(nz,4,7),
     $ saci(nz,4,7),sacw(nz,4,7),raci(nz,4,7),tacr(nz,4,7),raut(nz,4,7),
     $ racw(nz,4,7),sfw(nz,4,7),sfi(nz,4,7),gacs(nz,4,7),gacw(nz,4,7),
     $ gaci(nz,4,7),gacr(nz,4,7),gwet(nz,4,7),gaut(nz,4,7),racs(nz,4,7),
     $ sacr(nz,4,7),gfr(nz,4,7),smlt(nz,4,7),gmlt(nz,4,7),sdep(nz,4,7),
     $ ssub(nz,4,7),gsub(nz,4,7),pern(nz,4,7),d3ri(nz,4,7),d3ir(nz,4,7),
     $ d2sr(nz,4,7),d2rs(nz,4,7),gdry(nz,4,7),coc(nz,4,7),coe(nz,4,7),
     $ smf0(nz,4,7),qc0(nz,4,7),qr0(nz,4,7),qi0(nz,4,7),qs0(nz,4,7),
     $ qg0(nz,4,7),sqc0(nz,4,7),sqr0(nz,4,7),sqi0(nz,4,7),sqs0(nz,4,7),
     $ sqg0(nz,4,7),erns(nz,4,7),wgrs(nz,4,7),qsws(nz,4,7),
     $ tb0(nz,4),qb0(nz,4)
      common/bsts/ thom,tdw,tmlt,saut,saci,sacw,raci,tacr,raut,
     $             racw,sfw,sfi,gacs,gacw,gaci,gacr,gwet,gaut,racs,
     $             sacr,gfr,smlt,gmlt,sdep,ssub,gsub,pern,d3ri,d3ir,
     $             d2sr,d2rs,gdry,coc,coe,smf0,qc0,qr0,qi0,qs0,
     $             qg0,sqc0,sqr0,sqi0,sqs0,sqg0,erns,wgrs,qsws,
     $             tb0,qb0

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real    q1(nz,4,7),q2(nz,4,7),sub(nz,4,7),otr(nz,4,7),pir(nz),
     $        fmlt(nz,4,7),evap(nz,4,7),cond(nz,4,7),dep(nz,4,7)
      real    d1(nz),d2(nz),d3(nz),d4(nz),d5(nz)

      integer iaminf,k,k1,m,l
      real    a00,a1,a11
      real    afc,alf,als,alv,asc,avc,bs1,bs2,bs3,bs4
      real    cp,dtday,dtmas,period,pird,pird1

      save

c
       m=4
       period=6.
        iaminf=360
       a00=60./20.
       a11=60.
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.004e7
       avc=alv/cp
       afc=alf/cp
       asc=als/cp
       dtday=1./period
       dtmas=1./period
      do 10 k=1,nz
  10   pir(k)=pi(k)/512.
      do 15 k=2,nz-1
        pird=dtday*pir(k)*am(k)/(rho(k)*dz)
        pird1=dtday*am(k)/(rho(k)*dz*512.)
       sft(k)=sft(k)*pi(k)*dtday
       sfq(k)=-sfq(k)
cc     sfq(k)=-avc*sfq(k)*dtday
        d1(k)=(stf1(k)-stf1(k+1))*pird*a11
        d2(k)=-avc*(sqf1(k)-sqf1(k+1))*pird1*a11
         d3(k)=tls(k)*pi(k)*dtday
cc15     d4(k)=avc*qls(k)*dtday
         d4(k)=qls(k)
  15     d5(k)=qa(k)-taa(k)
      do 100 l=1,4
      do 100 k=1,nz
        pird=dtday/512.
       cond(k,l,m)=avc*coc(k,l,m)*a00
       evap(k,l,m)=-avc*coe(k,l,m)*a00
      bs1=avc*(coc(k,l,m)-coe(k,l,m))*a00
       dep(k,l,m)=asc*(gaut(k,l,m)+sdep(k,l,m)+gsub(k,l,m))*a00
       sub(k,l,m)=-asc*(pern(k,l,m))*a00
      bs2=asc*(gaut(k,l,m)-pern(k,l,m)+sdep(k,l,m)+gsub(k,l,m)-
     1    gwet(k,l,m)-gfr(k,l,m))*a00
       fmlt(k,l,m)=-afc*(smlt(k,l,m)+gmlt(k,l,m))*a00
      bs3=afc*(smlt(k,l,m)+gmlt(k,l,m))*a00
       otr(k,l,m)=afc*(sacw(k,l,m)+tacr(k,l,m)+sfw(k,l,m)+tmlt(k,l,m)+
     1            gacw(k,l,m)+gacr(k,l,m)+sacr(k,l,m)-d2rs(k,l,m))*a11
     2           -asc*(gwet(k,l,m)+gfr(k,l,m))*a00
      bs4=afc*(sacw(k,l,m)+tacr(k,l,m)+sfw(k,l,m)+tmlt(k,l,m)+
     1    gacw(k,l,m)+gacr(k,l,m)+sacr(k,l,m)-d2rs(k,l,m))*a11
       q1(k,l,m)=(bs1+bs2+bs4-bs3)*pird
cc     q2(k,l,m)=(bs1+bs2)*pird
       q2(k,l,m)=(bs1/avc+bs2/asc)/512.
       smf0(k,l,m)=smf0(k,l,m)*dtmas
        cond(k,l,m)=cond(k,l,m)*pird
        evap(k,l,m)=evap(k,l,m)*pird
        dep(k,l,m)=dep(k,l,m)*pird
        sub(k,l,m)=sub(k,l,m)*pird
        fmlt(k,l,m)=fmlt(k,l,m)*pird
        otr(k,l,m)=otr(k,l,m)*pird
  100  continue
       write(6,601) iaminf,m
      do 200 k1=1,nz
       k=nz+1-k1
  200 write(6,602) k,q1(k,1,m),q1(k,2,m),q1(k,3,m),q2(k,1,m),q2(k,2,m),
     1  q2(k,3,m),smf0(k,1,m),smf0(k,2,m),smf0(k,3,m)
       write(6,601) iaminf,m
      do 250 k1=1,nz
       k=nz+1-k1
  250 write(6,602) k,cond(k,1,m),cond(k,2,m),cond(k,3,m),evap(k,1,m),
     1  evap(k,2,m),evap(k,3,m),fmlt(k,1,m),fmlt(k,2,m),fmlt(k,3,m)
       write(6,601) iaminf,m
      do 270 k1=1,nz
       k=nz+1-k1
  270 write(6,602) k,dep(k,1,m),dep(k,2,m),dep(k,3,m),sub(k,1,m),
     1  sub(k,2,m),sub(k,3,m),otr(k,1,m),otr(k,2,m),otr(k,3,m)
      do 275 k1=1,nz
       k=nz+1-k1
       a1=-q2(k,1,4)-sfq(k)+d4(k)
  275 write(6,602) k,sfq(k),d4(k),a1,d5(k),d1(k),d2(k),d3(k),sft(k),
     1  sv(k)
  601 format(2x,'lev',4x,'q1(l=1)',4x,'q1(l=2)',4x,'q1(l=3)',4x,
     1 'q2(l=1)',4x,'q2(l=2)',4x,'q2(l=3)',4x,'mass(1)',4x,
     2 'mass(2)',4x,'mass(3)',10x,i4,'min',3x,i2,/)
  602 format(1x,i4,1x,9(1x,e10.3))
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rinit (irs)
c     ******   initialize all variables to zeros   ******
      parameter (NX=514,NZ=43,NT=2880,lay=88,ITT=244,nnt=481)
      parameter (nxi=nx-2,nxm=nx-1,nzm=nz-1)
      parameter (nx3=3*nx,nx16=16*nx,nz2=2*nz,nz4=4*nz,nz9=9*nz)
      parameter (nz13=13*nz,nz15=15*nz,nz20=20*nz,nz23=23*nz,nz28=28*nz)
      parameter (nb=nt*nxi,nb1=5*3*nt,nb8=12*nz+5*nt)
      parameter (nb5=178,nb6=193,nb7=2*nz+6,nb9=48*nz*28+2*nz*4)
      parameter (nb12=2*nxm*nzm,nb13=13*nz*28,nb14=3*nz+1)
      parameter (nb15=2*nz*28)
c      parameter (nb16=5*nx*nz*4)
      parameter (nz40=15*nz+5*nz*itt,NB24=NZ*21*14)
      PARAMETER (NB25=21*NZ*14+2*NZ,NB26=16*NZ*14+5*NZ*14)
      parameter (nx2=2*nx,nz17=17*nz,nnt3=3*nnt,nnt7=7*nnt)
      parameter (nb27=2*nz*4*7,NZ41=2*ITT*NZ+2*NZ,nz6=6*nz)
      parameter (nz42=2*ITT*nz+3*nz,NZ43=4*ITT*NX,NZ44=4*NX+2*NZ)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/mbudget/ acoc(nz17)
      common/tqave/ tavet(nz2)
      common/tbudget/ avett(nnt3)
      common/bsfc/ tsfc2(nx2)
      common/tbudget1/ t_sfcq(nnt7)
      COMMON/BSTS20/ STS20(NB27)

      COMMON/Q1Q2Z/ Q1Q2(NZ41)
      COMMON/Q1Q2T/ AQ1ZT(NZ4)
      COMMON/UPPER1/ T_ADJUST(NZ42)
      COMMON/BLS99/ THS(NZ43)
      COMMON/DINRAD/ P00(NZ44)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real    tcon(10)
      common/cont/ tcon
      real    tervr(9)
      common/rterv/ tervr
      real    trsnw(nb6)
      common/rsnw/ trsnw
      real    tb3cs(12)
      common/b3cs/ tb3cs
      real    rfa11(nb14)
      common/damp/ rfa11
      integer ivx(3),ivz(3)
      common/bxy/ ivx,ivz
      real    vx(8),vz(8)
      common/bx/ vx
      common/bz/ vz
      integer iib(3)
      common/bstart/ iib
      real    rrb(3),ao4(nb7),vb(14)
      common/rbstart/ rrb
      common/o4/ ao4
      common/bb/ vb
      real    v11(nx,nz),v12(nx,nz),v13(nx,nz),v14(nx,nz),
     $        v21(nx,nz),v22(nx,nz),v23(nx,nz),v24(nx,nz)
      common/b1t/ v11
      common/b1q/ v12
      common/b1c/ v13
      common/b1r/ v14
      common/b1u/ v21
      common/b1v/ v22
      common/b1w/ v23
      common/b1a/ v24
      real    vb1(nx,nz),vb2(nx,nz),vb3(nx,nz),vb4(nx,nz),
     $        vb5(nx,nz),vb6(nx,nz),vb7(nx,nz),vb8(nx,nz)
      common/b2t/ vb1
      common/b2q/ vb2
      common/b2c/ vb3
      common/b2r/ vb4
      common/b2u/ vb5
      common/b2v/ vb6
      common/b2w/ vb7
      common/b2a/ vb8
      real    qi (nx,nz),qs (nx,nz),qg (nx,nz),
     $        qi1(nx,nz),qs1(nx,nz),qg1(nx,nz)
      common/b1i/ qi
      common/b1s/ qs
      common/b1g/ qg
      common/b2i/ qi1
      common/b2s/ qs1
      common/b2g/ qg1
      real    v31(nx,nz),v32(nx,nz),v33(nx,nz),v34(nx,nz)
      common/bsat/ v31
      common/bsat1/ v32
      common/badv/ v33
      common/badv1/ v34
      real    sw1(nx,nz),sw2(nx,nz)
      common/slwave/ sw1,sw2
      real    vw(nx,nz),df0(nx,nz),trah(nx,nz),trav(nx,nz)
      common/bw/ vw,df0,trah,trav
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real    rst(nb),pcltop(nb),pclbot(nb),vw2(nb),vw3(nb1),VW22(NB)
c      real vw3(nb1)
      common/bw1/ rst,pcltop,pclbot 
      common/bw2/ vw2,VW22
      common/bw3/ vw3

      real    y4(nz9),y5(nz23),y6(nz13),x6(nx3),y8(nz15),y7(nz20),
     $        y9(nz20),ym(nz20),yb(nz40),ybu(nz28)
      common/b4/ y4
      common/b5/ y5
      common/b6/ y6,x6
      common/b8/ y8,y7
      common/b9/ y9,ym
      common/bb6/ yb
      common/bbtu/ ybu

      real    xa(nx16)
      common/ba/ xa

      real    s1(nb24),s2(nb25),s22(nb26)
      common/btt/ s1
      common/bcs/ s2
      common/bcss/ s22

      real    vrh1(nb8),vterv(6),vsnw(nb5)
      common/brh1/ vrh1
      common/bterv/ vterv
      common/bsnw/ vsnw

      real    sts3(nb9),sts4(nb13),stls1(nb15),rfft(nb12)
c      real sts5(nb16)
      common/bsts/ sts3
      common/bsts1/ sts4
c      common/bstsi/ sts5
      common/stls/ stls1
      common/rfft/ rfft

      real    sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)
      common/srflx/sfir,sfsw,shir,shsw,salpha,si
c
      COMMON/SUE/ PPRESS(NX,NZ)
      COMMON/SUE1/ TMODO(NX,NZ),QMODO(NX,NZ),RAINCO(NX,1)
      COMMON/SUE3/ T_BM(NB27)
      common/b66b/ s_dep(nz6)

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

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer i,k,n

      save

ccccccccccccccccccc
      DO I=1,NX
        RAINCO(I,1)=0.
        DO K=1,NZ
          ppress(i,k)=0.
          tmodo(i,k)=0.
          qmodo(i,k)=0.
        ENDDO
      ENDDO
      DO I=1,NZ41
        Q1Q2(I)=0.
      ENDDO
      DO I=1,NZ4
        AQ1ZT(I)=0.
      ENDDO

      DO I=1,NZ42
        T_ADJUST(I)=0.
      ENDDO
      DO I=1,NZ43
        THS(I)=0.
      ENDDO
      DO I=1,NZ44
        P00(I)=0.
      ENDDO
ccccccccccccccccccc
       do 2 i=1,nb8
    2   vrh1(i)=0.0
       do 3 k=1,nz20
        y7(k)=0.
        y9(k)=0.
    3   ym(k)=0.
       do 4 k=1,nz15
    4   y8(k)=0.
      DO 74 I=1,NB24
  74    S1(I)=0.
      DO 76 I=1,NB25
  76    S2(I)=0.
      DO 78 I=1,NB26
  78    S22(I)=0.
       do 7 n=1,nb9
    7   sts3(n)=0.
       do 8 n=1,nb13
    8   sts4(n)=0.
       do n=1,nb27
        T_BM(n)=0.
        STS20(n)=0.
       enddo
       do 9 n=1,nb15
    9   stls1(n)=0.

      do k=1,nz17
        acoc(k)=0.
      enddo
      do k=1,nz2
        tavet(k)=0.
      enddo
      do ni=1,nnt3
        avett(ni)=0.
      enddo
      do ni=1,nnt7
        t_sfcq(ni)=0.
      enddo
      do i=1,nx2
        tsfc2(i)=0.
      enddo
c
      if(irs.eq.1) return
       do 11 i=1,nx3
   11   x6(i)=0.0
       do 12 i=1,3
        ivx(i)=0
        ivz(i)=0
        iib(i)=0
   12   rrb(i)=0.
       do 13 i=1,8
        vx(i)=0.
   13   vz(i)=0.
       do 14 i=1,14
   14   vb(i)=0.
       do 15 k=1,nz
       do 15 i=1,nx
        v11(i,k)=0.
        v12(i,k)=0.
        v13(i,k)=0.
        v14(i,k)=0.
        v21(i,k)=0.
        v22(i,k)=0.
        v23(i,k)=0.
        v24(i,k)=0.
        vb1(i,k)=0.
        vb2(i,k)=0.
        vb3(i,k)=0.
        vb4(i,k)=0.
        vb5(i,k)=0.
        vb6(i,k)=0.
        vb7(i,k)=0.
        vb8(i,k)=0.
        v31(i,k)=0.
        v32(i,k)=0.
        v33(i,k)=0.
        v34(i,k)=0.
        qi(i,k)=0.
        qg(i,k)=0.
        qs(i,k)=0.
        qi1(i,k)=0.
        qg1(i,k)=0.
        qs1(i,k)=0.
        sw1(i,k)=0.
        sw2(i,k)=0.

c        q1_g_h(i,k)=0.
c        q1_g_v(i,k)=0.
c        q1_d_h(i,k)=0.
c        q1_d_v(i,k)=0.
c        q2_g_h(i,k)=0.
c        q2_g_v(i,k)=0.
c        q2_d_h(i,k)=0.
c        q2_d_v(i,k)=0.
c        q1_hyd(i,k)=0.
c        q2_hyd(i,k)=0.
c        q1_rad(i,k)=0.
c        q1a_g_h(i,k)=0.
c        q1a_g_v(i,k)=0.
c        q1a_d_h(i,k)=0.
c        q1a_d_v(i,k)=0.
c        q2a_g_h(i,k)=0.
c        q2a_g_v(i,k)=0.
c        q2a_d_h(i,k)=0.
c        q2a_d_v(i,k)=0.
c        q1a_hyd(i,k)=0.
c        q2a_hyd(i,k)=0.
c        q1a_rad(i,k)=0.

   15   vw(i,k)=0.
       do 16 k=1,nz9
   16   y4(k)=0.
       do 17 k=1,nz23
   17   y5(k)=0.
       do 18 k=1,nz13
   18   y6(k)=0.
       do 19 i=1,nx16
   19   xa(i)=0.
       do 20 k=1,nz40
   20   yb(k)=0.
       do 22 n=1,nb
        rst(n)=0.
        pcltop(n)=0. 
        pclbot(n)=0.
        VW22(N)=0.
   22   vw2(n)=0.
       do 23 n=1,nb1
   23   vw3(n)=0.
       do 24 n=1,6
   24   vterv(n)=0.
       do 25 n=1,nb5
   25   vsnw(n)=0.
       do 26 n=1,nb7
   26   ao4(n)=0.
       do 27 n=1,nb12
   27   rfft(n)=0.
       do 28 n=1,10
   28   tcon(n)=0.
       do 29 n=1,9
   29   tervr(n)=0.
       do 30 n=1,nb6
   30   trsnw(n)=0.
       do 31 n=1,12
   31   tb3cs(n)=0.
       do 32 n=1,nb14
   32   rfa11(n)=0.
       do 33 i=1,4
       salpha(i)=0.
       si(i)=0.
       do 33 n=1,lay
        sfsw(n,i)=0.
        sfir(n,i)=0.
        shsw(n,i)=0.
   33   shir(n,i)=0.
      do k=1,nz6
         s_dep(k)=0.
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine terp1 (n,x,f,w,y,int,tab,itab)
      dimension x(n),f(n),w(n),tab(3),itab(3)
      save
      call search (y,x,n,i)
      call interp (n,x,f,w,y,i,int,tab,itab)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine search (xbar,x,n,i)
c-----------------------------------------------------------------------
c--- Search the level i at which the interpolation begins            ---
c--- x increases with decreasing height                              ---
c-----------------------------------------------------------------------
      real    x(n)
      data b/.69314718/
      save
      if(xbar.gt.x(2)) go to 101
      i=1
      return
  101 continue
      if(xbar.lt.x(n-1)) go to 102
      i=n-1
      return
  102 continue
      m=int((log(float(n)))/b)
      i=2**m
      if(i.ge.n) i=i/2
      k=i
      nm1=n-1
  103 continue
      k=k/2
      if(xbar.ge.x(i)) go to 104
      i=i-k
      go to 103
  104 continue
      if(xbar.le.x(i+1)) return
      i=min0(i+k,nm1)
      go to 103
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interp (n,x,f,w,y,i,int,tab,itab)
      dimension x(n),f(n),w(n),tab(3),itab(3)
      save
      ii(i)=(i-1)*int+1
      flk=x(i+1)-x(i)
      flp=x(i+1)-y
      fl0=y-x(i)
      i0=ii(i)
      ip=i0+int
      if(itab(1).ne.1) go to 102
      a=(w(i0)*flp**3+w(ip)*fl0**3)/(6.*flk)
      b=(f(ip)/flk-w(ip)*flk/6.)*fl0
      c=(f(i0)/flk-w(i0)*flk/6.)*flp
      tab(1)=a+b+c
  102 if(itab(2).ne.1) go to 104
      a=(w(ip)*fl0**2-w(i0)*flp**2)/(2.*flk)
      b=(f(ip)-f(i0))/flk
      c=(w(i0)-w(ip))*flk/6.
      tab(2)=a+b+c
  104 if(itab(3).ne.1) go to 106
      tab(3)=(w(i0)*flp+w(ip)*fl0)/flk
  106 return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coeff (n,x,f,w,iop,int,wk)
c-----------------------------------------------------------------------
c--- Compute the coefficients for the cubic spline interpolation     ---
c-----------------------------------------------------------------------
c--- Inputs:  x --- pressure at the model level                      ---
c---          f --- variables being interpolated                     ---
c--- Outputs: w --- coefficients                                     ---
c---          wk--- stores couple of coefficients at each point      ---
c-----------------------------------------------------------------------
      dimension x(n),f(n),w(n),wk(n,4),iop(2)
      save
      ii(i)=(i-1)*int+1
      j0=1
      do i=2,n
        jm=j0
        j0=j0+int
        wk(i,1)=x(i)-x(i-1)
        wk(i,2)=(f(j0)-f(jm))/wk(i,1)
        wk(i,3)=wk(i,1)/6.
        wk(i,1)=wk(i,1)/3.
      end do
      nn=n
      mk=iop(1)
      ml=iop(2)
      go to (102,103,104,105) ,mk
  102 continue
      wk(2,2)=wk(3,2)-wk(2,2)-wk(2,3)*w(1)
      wk(2,3)=0.
      wk(2,1)=wk(2,1)+wk(3,1)
      i1=2
      nn=nn-1
      go to 106
  103 continue
      wk(1,2)=wk(2,2)-w(1)
      wk(2,2)=wk(3,2)-wk(2,2)
      wk(1,3)=0.
      wk(1,1)=wk(2,1)
      wk(2,1)=wk(2,1)+wk(3,1)
      i1=1
      go to 106
  104 continue
      y2=wk(2,2)
      b2=wk(2,1)
      wk(2,2)=wk(3,2)-wk(2,2)
      wk(2,1)=wk(3,1)+wk(2,1)
      i1=2
      nn=nn-1
      go to 106
  105 continue
      a12=x(1)-x(2)
      a13=x(1)-x(3)
      a14=x(1)-x(4)
      a23=x(2)-x(3)
      a24=x(2)-x(4)
      a34=x(3)-x(4)
      j1=1
      j2=j1+int
      j3=j2+int
      j4=j3+int
      w(1)=(1./a12+1./a13+1./a14)*f(j1)-
     1     a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-
     2     a12*a13/(a14*a24*a34)*f(j4)
      go to 103
  106 continue
      i2=n-2
      do i=3,i2
        wk(i,2)=wk(i+1,2)-wk(i,2)
        wk(i,1)=wk(i+1,1)+wk(i,1)
      end do
      in=ii(n)
      go to (108,109,110,111) ,ml
  108 continue
      wk(n-1,2)=wk(n,2)-wk(n-1,2)-wk(n,3)*w(in)
      wk(n,3)=0.
      wk(n-1,1)=wk(n-1,1)+wk(n,1)
      nn=nn-1
      go to 112

  109 continue
      wk(n-1,2)=wk(n,2)-wk(n-1,2)
      wk(n,2)=-wk(n,2)+w(in)
      wk(n-1,1)=wk(n-1,1)+wk(n,1)
      wk(1,4)=0.
      go to 112

  110 continue
      wk(n-1,2)=wk(n,2)-wk(n-1,2)
      wk(n,2)=y2-wk(n,2)
      wk(n-1,1)=wk(n-1,1)+wk(n,1)
      wk(n,1)=wk(n,1)+b2
      wk(1,4)=wk(2,3)
      go to 112

  111 continue
      a12=x(n)-x(n-1)
      a13=x(n)-x(n-2)
      a14=x(n)-x(n-3)
      a23=x(n-1)-x(n-2)
      a24=x(n-1)-x(n-3)
      a34=x(n-2)-x(n-3)
      j1=in
      j2=j1-int
      j3=j2-int
      j4=j3-int
      w(in)=(1./a12+1./a13+1./a14)*f(j1)-
     1      a13*a14/(a12*a23*a24)*f(j2)+a12*a14/(a13*a23*a34)*f(j3)-
     2      a12*a13/(a14*a24*a34)*f(j4)
      go to 109

  112 continue
      ii1=ii(i1)
      call trip (nn,wk(i1,3),wk(i1,1),wk(i1+1,3),wk(i1,2),w(ii1),int)
      go to (114,114,113,114) ,mk

  113 continue
      w(1)=w(in)
  114 continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trip (n,a,b,c,y,z,int)
      real    a(n),b(n),c(n),y(n),z(n)
      save
      ii(i)=(i-1)*int+1
      bn=b(n)
      yn=y(n)
      v=c(n)
      y(1)=y(1)/b(1)
      a(1)=a(1)/b(1)
      b(1)=c(1)/b(1)
      nm2=n-2
      do j=2,nm2
        den=b(j)-a(j)*b(j-1)
        b(j)=c(j)/den
        y(j)=(y(j)-a(j)*y(j-1))/den
        a(j)=-a(j)*a(j-1)/den
        bn=bn-v*a(j-1)
        yn=yn-v*y(j-1)
        v=-v*b(j-1)
      end do
      den=b(n-1)-a(n-1)*b(n-2)
      b(n-1)=(c(n-1)-a(n-1)*a(n-2))/den
      y(n-1)=(y(n-1)-a(n-1)*y(n-2))/den
      bn=bn-v*a(n-2)
      yn=yn-v*y(n-2)
      v=a(n)-v*b(n-2)
      nm1=n-1
      in=ii(n)
      inm=ii(nm1)
      z(in)=(yn-v*y(nm1))/(bn-v*b(nm1))
      z(inm)=y(nm1)-b(nm1)*z(in)
      do j=2,nm1
        k=n-j
        ik=ii(k)
        z(ik)=y(k)-b(k)*z(ik+int)-a(k)*z(in)
      end do
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rainstat (rst)

      implicit none
      integer nx,nt,nstr,nxi
      PARAMETER (NX=514,NT=2880,NSTR=100)
      PARAMETER (NXI=NX-2)

      real    rst(nt,nxi)

      real    area(nstr),accum(nstr),xinter(nstr)
      integer ik,it,ix
      real    ri,rmin,sum1,sum2

      save
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ----------------------------------------------------------------------
      RMIN = 0.01
      DO 1 IK = 1,NSTR
        XINTER(IK) = FLOAT(IK)
        AREA(IK)=0.
        ACCUM(IK)=0.
    1 CONTINUE
      DO 100 IT = 1,NT
      DO 100 IX = 1,NXI
        RI=RST(IT,IX)
        IF (RI .GE. RMIN) THEN
          IF (RI .LT. XINTER(1)) THEN
            ACCUM(1) = ACCUM(1) + RI
            AREA(1) = AREA(1) + 1.
          ENDIF 
          DO 10 IK = 2,NSTR
            IF (RI .GE. XINTER(IK-1) .AND. RI .LT. XINTER(IK)) THEN
              ACCUM(IK) = ACCUM(IK) + RI
              AREA(IK) = AREA(IK) + 1.
            ENDIF
   10     CONTINUE
        ENDIF
  100 CONTINUE
      SUM1 = 0.
      SUM2 = 0.
      DO IK = 1,NSTR
       SUM1 = SUM1 + ACCUM(IK)
       SUM2 = SUM2 + AREA(IK)
      ENDDO
      DO IK = 1,NSTR
        IF (SUM1 .EQ. 0.) THEN
          ACCUM(IK) = 0.
        ELSE
          ACCUM(IK) = 100.*ACCUM(IK)/SUM1
        ENDIF
        IF (SUM2 .EQ. 0.) THEN
          AREA(IK) = 0.
        ELSE
          AREA(IK) = 100.*AREA(IK)/SUM2
        ENDIF
      ENDDO
      DO 200 IX = 1,NSTR
        WRITE (6,101) IX,XINTER(IX),ACCUM(IX)
        WRITE (6,102) IX,AREA(IX)
  200 CONTINUE
  101 FORMAT(12HACCUM(IX) = ,I8,4X,2F12.3)
  102 FORMAT(12HAREA(IX)  = ,I8,4X,2F12.3)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TIMING(T)
      EXTERNAL ETIME         !sgi
      REAL*4 ETIME,TARRAY(2) !sgi
      T=ETIME(TARRAY)        !sgi
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine budget_m (dt_stat)
      parameter (NX=514,NZ=43)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/mbudget/ acoc(nz),acoe(nz),agaut(nz),asdep(nz),agsub(nz),
     1   apern(nz),aqls(nz),altrans(nz),actrans(nz),total_mq(nz),
     2   agwet(nz),afgr(nz),avelw(nz),avesw(nz),amelt(nz),aother(nz),
     3   total_mt(nz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
     $ tb0(nz,4),qb0(nz,4)
      common/stls/ srsw(nz,4,7),srlw(nz,4,7)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       a22=1./(nx-2)
       a11=dt_stat
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.004e7
       avc=alv/cp
       afc=alf/cp
       asc=als/cp
      do k=1,nz
        acoc(k)=acoc(k)+coc(k,1,4)*a22
        acoe(k)=acoe(k)+coe(k,1,4)*a22
        agaut(k)=agaut(k)+gaut(k,1,4)*a22
        asdep(k)=asdep(k)+sdep(k,1,4)*a22
        agsub(k)=agsub(k)+gsub(k,1,4)*a22
        apern(k)=apern(k)+pern(k,1,4)*a22
        agwet(k)=agwet(k)+gwet(k,1,4)*a22
        afgr(k)=afgr(k)+gfr(k,1,4)*a22
        amelt(k)=amelt(k)+(smlt(k,1,4)+gmlt(k,1,4))*a22
        aother(k)=aother(k)+(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     1                      tmlt(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     2                      sacr(k,1,4)-d2rs(k,1,4))*a11*a22
        total_mq(k)=total_mq(k)+(coc(k,1,4)-coe(k,1,4)
     1            +gaut(k,1,4)+sdep(k,1,4)+gsub(k,1,4)-pern(k,1,4)
     2            -gwet(k,1,4)-gfr(k,1,4))*a22
        total_mt(k)=total_mt(k)+avc*a22*(coc(k,1,4)-coe(k,1,4))
     1              +asc*a22*(gaut(k,1,4)+sdep(k,1,4)+gsub(k,1,4)
     2                    -pern(k,1,4)-gwet(k,1,4)-gfr(k,1,4))
     3              -afc*a22*(smlt(k,1,4)+gmlt(k,1,4))
     4              +afc*a22*(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     5                    tmlt(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     6                    sacr(k,1,4)-d2rs(k,1,4))*a11
        avelw(k)=avelw(k)+srlw(k,1,4)*a11*a22
        avesw(k)=avesw(k)+srsw(k,1,4)*a11*a22
      end do
      return
      end
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine budget_mt (dt_stat,total_mq,total_mt)
      parameter (NX=514,NZ=43)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension total_mq(nz),total_mt(nz)
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
     $ tb0(nz,4),qb0(nz,4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       a22=1./(nx-2)
       a11=dt_stat
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.004e7
       avc=alv/cp
       afc=alf/cp
       asc=als/cp
      do k=1,nz
        total_mq(k)=(coc(k,1,4)-coe(k,1,4)
     1            +gaut(k,1,4)+sdep(k,1,4)+gsub(k,1,4)-pern(k,1,4)
     2            -gwet(k,1,4)-gfr(k,1,4))*a22
        total_mt(k)=avc*a22*(coc(k,1,4)-coe(k,1,4))
     1              +asc*a22*(gaut(k,1,4)+sdep(k,1,4)+gsub(k,1,4)
     2                    -pern(k,1,4)-gwet(k,1,4)-gfr(k,1,4))
     3              -afc*a22*(smlt(k,1,4)+gmlt(k,1,4))
     4              +afc*a22*(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     5                    tmlt(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     6                    sacr(k,1,4)-d2rs(k,1,4))*a11
      end do
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine var_ave (t,tave)
      parameter (NZ=43)
c***********************************************************************
c*** this program computes the quantities displayed in fig.6 of      ***
c*** sui et al. (1994) jas v51. 711-728                              ***
c***-----------------------------------------------------------------***
c*** input variables:                                                ***
c***       kk........which level to display the thee field           ***
c***       t=qv.........density field                                ***
c***-----------------------------------------------------------------***
c*** output variables:                                               ***
c***       tave .....domain averaged and density weighted water vapor***
c***********************************************************************
      dimension t(nz)
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     1   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     2   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
c compute vertically averaged field
      tave  =0.0
      do 500 k=2,nz-1
        aa    =dz*rho(k)/am(k)
        tave  =tave  +aa*t(k)
 500  continue
      tave=10.*tave
      print*,'denq (mm)=',tave
      return
      end
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine var_avet (t,tave)
c      parameter (NX=514,NZ=43,ITT=244)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc*** input variables:                                                ***
cc***       t.........field                                           ***
cc***-----------------------------------------------------------------***
cc*** output variables:                                               ***
cc***       tave  ....domain averaged and density weighted temperature***
cc***********************************************************************
c      real t(nz),t1(nz)
c      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
c      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
c     1   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
c     2   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
c      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
c     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
c      save
cc  compute vertically averaged field
c      denave = 0.0
c      tave=0.0
c      do k=2,nz-1
c        aa=dz*rho(k)/am(k)
c        denave=denave+aa
c        tave  =tave  +aa*t1(k)
c      enddo
c      tave  =tave/denave
c        print *,'tave (k)= ',tave
c      return
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dom_ave (time,tave,denq,thee)
      parameter (NX=514,NZ=43)
c***********************************************************************
c*** This program computes the quantities displayed in Fig.6 of      ***
c*** Sui et al. (1994) JAS V51. 711-728                              ***
c***-----------------------------------------------------------------***
c*** Input variables:                                                ***
c***       kk........which level to display the thee field           ***
c***       t.........density field                                   ***
c***       the.......potential temperature field purtabation         ***
c***       qv........water vapor field purtabation                   ***
c***       p.........purtabation pressure field                      ***
c***       pi........base state pi field                             ***
c***       ta........base potential temperature                      ***
c***       qa........base q (kg/kg)                                  ***
c***       id........=1 when 2D model is used (y=1)                 ***
c***-----------------------------------------------------------------***
c*** Output variables:                                               ***
c***       tave  ....domain averaged and density weighted temperature***
c***       denq.....domain averaged and density weighted water vapor***
c***       thee......equivalent potential temperature                ***
c***********************************************************************
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      COMMON/BB/ DT,D2T,RIL2,F5,RD1,RD2,BOUND,AL,CP,RA,CK,CE,EPS,PSFC
      COMMON/B2T/ DPT1(NX,NZ)
      COMMON/B2Q/ DQV1(NX,NZ)
      COMMON/B5/ TB(NZ),QB(NZ),RHO1(NZ),RHO(NZ),BA(NZ),BB(NZ),TA(NZ),
     1   QA(NZ),TA1(NZ),QA1(NZ),COEF(NZ),C1(NZ),C2(NZ),C3(NZ),AM(NZ),
     2   AM1(NZ),UB(NZ),VB(NZ),WB(NZ),UB1(NZ),VB1(NZ),RRHO(NZ),RRHO1(NZ)
      COMMON/B6/ TLS(NZ),QLS(NZ),FD(NZ),FE(NZ),P0(NZ),PI(NZ),F0(NZ),
     1   ST(NZ),SV(NZ),SQ(NZ),SC(NZ),SE(NZ),SQA(NZ),RI(NX),AR(NX),RX(NX)
      common/tvertical/ denave
c local variables
      real    t(1000)
      save
      sx  = real(nx-2)
      do 100 k=2,nz-1
        t(k)=ta1(k)*PI(k)
 100  continue
c compute vertically averaged field
      denave=0.0
      tave  =0.0
      denq  =0.0
      dzt   =0.0
      do 500 k=2,nz-1
        aa    =dz*rho(k)/am(k)
        dzt   =dzt+dz/am(k)
        denave=denave+aa
        tave  =tave  +aa*t(k)
        denq  =denq  +aa*qa1(k)
 500  continue
      print *,'DZT=',DZT
      print *,'denave=',denave
      tave  =tave/denave
c     convert DENQ to mm
      denq  =10.0*denq
c compute THETAE at 950 mb
      kk=3
      thd=0.0
      do 600 i=2,nx-1
        t1=ta1(kk)+dpt1(i,kk)
        q1=qa1(kk)+dqv1(i,kk)
        thd=thd+t1*exp(1.0e3*2.4925263*q1/(t1*pi(kk)))
 600  continue
      thee=thd/sx
        print *,'time            = ',time
        print *,'temperature ave = ',tave
        print *,'water vapor ave = ',denq
        print *,'thret at 950 mb = ',thee
      return
      end

      subroutine consat (rho)
c    (r&h)  specify some constants in satice routine   ******
      parameter (NZ=43,NT=2880,nb=10*nz+5*nt)
      common/iceopt/ ice913,ilif
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
      common/rterv/ zrc,zgc,zsc,vrc0,vrc1,vrc2,vrc3,vgc,vsc
      common/rsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,rn2,bnd2,rn3,rn4,
     $  rn5,rn50,rn51,rn52,rn53,rn6,rn60,rn61,rn62,rn63,rn7,rn8,rn9,
     $  rn10,rn101,rn102,rn10a,rn10b,rn10c,rn11,rn12,rn12a(31),
     $  rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn171,rn172,rn17a,rn17b,
     $  rn17c,rn18,rn18a,rn19,rn191,rn192,rn19a,rn20,rn20a,rn20b,rn30,
     $  rn30a,rn21,bnd21,rn22,rn23,rn231,rn232,rn25,rn25a(31),rn31,beta,
     $  rn32,rn33,rn331,rn332,rn34,rn35
      common/brh1/ srr(nz),qrr(nz),z1(nb)
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/size/ tnw,tns,tng,roqs,roqg,roqr
      real    a1(31),a2(31),rho(nz)
      data a1/.7939e-7,.7841e-6,.3369e-5,.4336e-5,.5285e-5,.3728e-5,
     1   .1852e-5,.2991e-6,.4248e-6,.7434e-6,.1812e-5,.4394e-5,.9145e-5,
     2   .1725e-4,.3348e-4,.1725e-4,.9175e-5,.4412e-5,.2252e-5,.9115e-6,
     3   .4876e-6,.3473e-6,.4758e-6,.6306e-6,.8573e-6,.7868e-6,.7192e-6,
     4   .6513e-6,.5956e-6,.5333e-6,.4834e-6/
      data a2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
      save
c     ******************************************************************
      cpi=4.*atan(1.)
      cpi2=cpi*cpi
c      grvt=980.
      c76=7.66
      c358=35.86
      c172=17.26939
      c409=4098.026
      c218=21.87456
      c580=5807.695
      c610=6.1078e3
      c149=1.496286e-5
      c879=8.794142
      c141=1.4144354e7
c     ***************
        tca=2.43e3
        dwv=.226
        dva=1.718e-4
        amw=18.016
        ars=8.314e7
      t0=273.16
      t00=238.16
      alv=2.5e10
      alf=3.336e9
      als=2.8336e10
      avc=alv/cp
      afc=alf/cp
      asc=als/cp
      rw=4.615e6
      cw=4.187e7
      ci=2.093e7
c***   define the density and size distribution of precipitation
      roqr=1.
      tnw=.08
        roqs=.1
        tns=.16
cfred   tns=1.
          roqg=.4
cfred          tng=.04
        tng=.08
c***   define the coefficients used in terminal velocity
       ag=351.2
       bg=.37
          as=78.63154
          bs=.11
c         AS=152.93
c         BS=.25
         aw=2115.
         bw=.8
c
       if (ice913 .eq. 1) then
          t00=238.16
          tns=0.08
          tng=.04
          ag=372.3
          as=78.63154
          bs=.11
       endif
c
       bgh=.5*bg
       bsh=.5*bs
       bwh=.5*bw
       bgq=.25*bg
       bsq=.25*bs
       bwq=.25*bw
      ga3=2.
      ga4=6.
      ga5=24.
      ga6=120.
      ga7=720.
      ga8=5040.
      ga9=40320.
        ga4g=11.63177
        ga3g=3.3233625
        ga5gh=1.608355
        if(bg.eq.0.37) ga4g=9.730877
        if(bg.eq.0.37) ga3g=2.8875
        if(bg.eq.0.37) ga5gh=1.526425
          ga3d=2.54925
          ga4d=8.285063
          ga5dh=1.456943
          if(bs.eq.0.57) ga3d=3.59304
          if(bs.eq.0.57) ga4d=12.82715
          if(bs.eq.0.57) ga5dh=1.655588
          if(bs.eq.0.11) ga3d=2.218906
          if(bs.eq.0.11) ga4d=6.900796
          if(bs.eq.0.11) ga5dh=1.382792

cccccc             rutledge and hobbs, 1984   cccccccccccccccccccccccccc
          ga6d=144.93124
        ac1=as
        ac2=ag
       zrc=(cpi*roqr*tnw)**0.25
       zsc=(cpi*roqs*tns)**0.25
       zgc=(cpi*roqg*tng)**0.25
       vrc0=-26.7
       vrc1=20600./zrc
       vrc2=-204500./(zrc*zrc)
       vrc3=906000./(zrc*zrc*zrc)
       vsc=ac1*ga4d/(6.*zsc**bs)
       vgc=ac2*ga4g/(6.*zgc**bg)
cs      cd1=6.e-1
cs      cd2=4.*grvt/(3.*cd1)
cs     vgc=ga4g*sqrt(cd2*roqg/zgc)/6.
c     ****************************
      rn1=9.4e-15
      rn2=1.e-3
       bnd2=2.0e-3
c
       if (ice913 .eq. 1)  bnd2=1.5e-3
c
       esi=.1
      rn3=.25*cpi*tns*ac1*esi*ga3d
       esc=1.
      rn4=.25*cpi*esc*tns*ac1*ga3d
       eri=1.
CLIN   ERI=1.
       ERI=0.1
c
       if (ice913 .eq. 1)  ERI=1.
c
      rn5=.25*cpi*eri*tnw
       rn50=-.267e2*ga3
       rn51=5.15e3*ga4
       rn52=-1.0225e4*ga5
       rn53=7.55e3*ga6
       ami=1./(24.*6.e-9)
      rn6=cpi2*eri*tnw*roqr*ami
       rn60=-.267e2*ga6
       rn61=5.15e3*ga7
       rn62=-1.0225e4*ga8
       rn63=7.55e3*ga9
       esr=.5
c
       if (ice913 .eq. 1)  esr=1.
c
      rn7=cpi2*esr*tnw*tns*roqs
       esr=1.
      rn8=cpi2*esr*tnw*tns*roqr
       egs=.1
      rn9=cpi2*egs*tns*tng*roqs
      rn10=4.*tns
       rn101=.65
       rn102=.44*sqrt(ac1/dva)*ga5dh
       rn10a=alv*als*amw/(tca*ars)
       rn10b=alv/tca
       rn10c=ars/(dwv*amw)
      rn11=2.*cpi*tns*tca/alf
c HFO  
c	ami50=4.8e-7
c KFL (1993)
       ami50=4.8e-7*(100./50.)**3
       ami40=2.46e-7
       ami40=2.46e-7*.5**3
c       ami40=3.84e-9
       eiw=1.
       ui50=100.
c HFO
c	ri50=5.e-3
c KFL (1993)
       ri50=2.*5.e-3
       cmn=1.05e-15
      rn12=cpi*eiw*ui50*ri50*ri50
      do 10 k=1,31
        y1=1.-a2(k)
       rn13(k)=a1(k)*y1/(ami50**y1-ami40**y1)
       rn12a(k)=rn13(k)/ami50
       rn12b(k)=a1(k)*ami50**a2(k)
       rn25a(k)=a1(k)*cmn**a2(k)
   10 continue
       egc=1.
      rn14=.25*cpi*egc*ac2*tng*ga3g
       egi=.1
      rn15=.25*cpi*egi*tng*ac2*ga3g
       egi=1.
      rn15a=.25*cpi*egi*tng*ac2*ga3g
       egr=1.
      rn16=cpi2*egr*tng*tnw*roqr
      rn171=2.*cpi*tng*alv*dwv
       rn172=2.*cpi*tng*tca
       rn17a=.31*ga5gh*sqrt(ac2/dva)
       rn17b=cw-ci
       rn17c=cw
       apri=.66
       bpri=1.e-4
CTAO
       bpri=.5*bpri
c
       if (ice913 .eq. 1)  bpri=1.e-4
c
      rn18=20.*cpi2*bpri*tnw*roqr
       rn18a=apri
      rn19=2.*cpi*tng*tca/alf
       rn191=.78
       rn192=.31*ga5gh*sqrt(ac2/dva)
       rn19a=cw/alf
      rn20=2.*cpi*tng
       rn20a=als*als*amw/(tca*ars)
       rn20b=als/tca
      rn30=2.*cpi*tng
       rn30a=alv*alv*amw/(tca*ars)
      rn21=1.e-3
       bnd21=1.5e-3
       erc=1.
      rn22=.25*cpi*erc*tnw
      rn23=2.*cpi*tnw
       rn231=.78
       rn232=.31*ga3*sqrt(3.e3/dva)
       cn0=1.e-8

c      rn25=cn0	!tao-new
       rn25=cn0

      rn31=1.e-17
       beta=-.6
      rn32=4.*51.545e-4
      rn33=4.*tns
       rn331=.65
       rn332=.44*sqrt(ac1/dva)*ga5dh
CTAO
       esc=0.5
c
       if (ice913 .eq. 1)  esc=1.
c
       amc=1./(24.*4.e-9)
      rn34=cpi2*esc*amc*ac1*roqs*tns*ga6d
      rn35=alv*alv/(cp*rw)
c     ****************************
      do 20 k=1,nz
       srr(k)=1./sqrt(rho(k))
   20  qrr(k)=sqrt(srr(k))
      return
      end
      subroutine consatl (rho)
C     (LIN) SPECIFY SOME CONSTANTS IN SATICE ROUTINE   ******
      parameter (NZ=43,NT=2880)
      parameter (nb=10*nz+5*nt)

      real    rho(nz)
      real    srr(nz),qrr(nz),z1(nb)
      common/brh1/ srr,qrr,z1

      common/ilin/ lin

      common/iceopt/ ice913,ilif

      real    dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc

      real    c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
      common/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141

      real    ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
      common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq

      real    tnw,tns,tng,roqs,roqg,roqr
      common/size/ tnw,tns,tng,roqs,roqg,roqr

      COMMON/BTERV/ ZRC,ZGC,ZSC,VRC,VGC,VSC
      COMMON/BSNW/ ALV,ALF,ALS,T0,T00,AVC,AFC,ASC,RN1,BND1,RN2,BND2,
     1   RN3,RN4,RN5,RN6,RN7,RN8,RN9,RN10,RN101,RN10A,RN11,RN11A,
     2   RN12,RN12A(31),RN12B(31),RN13(31),RN14,RN15,RN15A,RN16,RN17,
     3   RN17A,RN17B,RN17C,RN18,RN18A,RN19,RN19A,RN19B,RN20,RN20A,RN20B,
     4   BND3,RN21,RN22,RN23,RN23A,RN23B,RN25,RN25A(31),RN30A,RN30B,
     5   RN30C,RN31,BETA,RN32
      DIMENSION A1(31),A2(31)
      DATA A1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,
     1   .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,
     2   .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,
     3   .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,
     4   .6513E-6,.5956E-6,.5333E-6,.4834E-6/
      DATA A2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,
     1   .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,
     2   .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,
     3   .4382,.4361/
C     ******************************************************************
        CPI=4.*ATAN(1.)
        CPI2=CPI*CPI
        GRVT=980.
       TCA=2.43E3
       DWV=.226
       AMW=18.016
       ARS=8.314E7
      T0=273.16
      T00=238.16
      ALV=2.5E10
      ALF=3.336E9
      ALS=2.8336E10
       AVC=ALV/CP
       AFC=ALF/CP
       ASC=ALS/CP
      RW=4.615E6
      CW=4.187E7
      CI=2.093E7
      C76=7.66
      C358=35.86
      C172=17.26939
      C409=4098.026
      C218=21.87456
      C580=5807.695
      C610=6.1078E3
      C149=1.496286E-5
      C879=8.794142
      C141=1.4144354E7
C***   DEFINE THE COEFFICIENTS USED IN TERMINAL VELOCITY
       IF (LIN .EQ. 1) THEN
        AG=1400.
        BG=.5
       ELSE
        AG=351.2
        BG=.37
       ENDIF
CRH    AG=351.2
CRH    BG=.37
        AS=152.93
        BS=.25
CYOSH   AS=4.0858E2
CYOSH   BS=.57
CRH     AS=68.63154
CRH     BS=.11
         AW=2115.
         BW=.8
c
       if (ice913 .eq. 1 .AND. LIN .EQ. 0) then
          ag=372.3
          bs=.11
       endif
c
       BGH=.5*BG
       BSH=.5*BS
       BWH=.5*BW
       BGQ=.25*BG
       BSQ=.25*BS
       BWQ=.25*BW
C***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
      ROQR=1.
      TNW=.08
CLORD TNW=.22
        ROQS=.084
CYOSH   ROQS=.084
        TNS=.0765
CRH     TNS=.04
CYOSH   TNS=.0765
CFRED   TNS=1.
c
CRH       ROQG=.4
CYOSH     ROQG=.9
CLORD     ROQG=.3
c
CRH       TNG=.04
CFRED     TNG=.22
         IF (LIN .EQ. 1) THEN
          ROQG=.9
          TNG=.002
         ELSE
          ROQG=.4
          TNG=.1
         ENDIF
      GA3B=4.6941552
      GA4B=17.83779
      GA6B=496.6041
      GA5BH=1.827363
        GA4G=11.63177
        GA3G=3.3233625
        GA5GH=1.608355
        IF(BG.EQ.0.37) GA4G=9.730877
        IF(BG.EQ.0.37) GA3G=2.8875
        IF(BG.EQ.0.37) GA5GH=1.526425
          GA3D=2.54925
          GA4D=8.285063
          GA5DH=1.456943
          IF(BS.EQ.0.57) GA3D=3.59304
          IF(BS.EQ.0.57) GA4D=12.82715
          IF(BS.EQ.0.57) GA5DH=1.655588
          IF(BS.EQ.0.11) GA3D=2.218906
          IF(BS.EQ.0.11) GA4D=6.900796
          IF(BS.EQ.0.11) GA5DH=1.382792
CCCCCC        LIN ET AL., 1983 OR LORD ET AL., 1984   CCCCCCCCCCCCCCCCC
      AC1=AW
      CC1=AS
      CD1=6.E-1
       CD2=4.*GRVT/(3.*CD1)
      ZRC=(CPI*ROQR*TNW)**0.25
      ZSC=(CPI*ROQS*TNS)**0.25
      ZGC=(CPI*ROQG*TNG)**0.25
      VRC=AC1*GA4B/(6.*ZRC**BW)
      VSC=CC1*GA4D/(6.*ZSC**BS)
      IF (LIN .EQ. 1) THEN
         VGC=GA4G*SQRT(CD2*ROQG/ZGC)/6.
      ELSE
         VGC=AG*GA4G/(6.*ZGC**BG)
      ENDIF
CRH   VGC=AG*GA4G/(6.*ZGC**BG)
CFRED VGC=AG*GA4G/(6.*ZGC**BG)
C     ****************************
      RN1=1.E-3
      RN2=1.E-3
      IF (LIN .EQ. 1) THEN
       BND1=7.5E-4
       BND2=1.25E-3
      ELSE
       BND1=1.E-3
       BND2=2.00E-3
c
       if (ice913 .eq. 1)  bnd2=1.5e-3
c
      ENDIF
CYOSH  BND2=3.E-3

      RN3=.25*CPI*TNS*CC1*GA3D
       ESW=1.
      RN4=.25*CPI*ESW*TNS*CC1*GA3D
CLIN       ERI=1.
       ERI=0.1
      RN5=.25*CPI*ERI*TNW*AC1*GA3B
       AMI=1./(24.*4.19E-10)
      RN6=CPI2*ERI*TNW*AC1*ROQR*GA6B*AMI
       ESR=1.
ctao
       ESR=0.5*ESR
      RN7=CPI2*ESR*TNW*TNS*ROQS
       ESR=1.
      RN8=CPI2*ESR*TNW*TNS*ROQR
       ESR=0.1*ESR
      RN9=CPI2*TNS*TNG*ROQS
      RN10=2.*CPI*TNS
       RN101=.31*GA5DH*SQRT(CC1)
       RN10A=ALS*ALS/RW
      RN11=2.*CPI*TNS/ALF
       RN11A=CW/ALF
CLIN   AMI50=4.8E-7
       AMI50=4.8E-7*(100./50.)**3
       AMI40=3.84E-9
       EIW=1.
       UI50=100.
CLIN   RI50=5.E-3
       RI50=2.*5.E-3
       CMN=1.05E-15
      RN12=CPI*EIW*UI50*RI50**2
      DO 10 K=1,31
       Y1=1.-A2(K)
      RN13(K)=A1(K)*Y1/(AMI50**Y1-AMI40**Y1)
       RN12A(K)=RN13(K)/AMI50
       RN12B(K)=A1(K)*AMI50**A2(K)
       RN25A(K)=A1(K)*CMN**A2(K)
   10 CONTINUE
       EGW=1.
      RN14=.25*CPI*EGW*TNG*GA3G*SQRT(CD2*ROQG)
CLIN    EGI=.1
        EGI=0.05
      RN15=.25*CPI*EGI*TNG*GA3G*SQRT(CD2*ROQG)
CLIN   EGI=1.
       EGI=0.5
      RN15A=.25*CPI*EGI*TNG*GA3G*SQRT(CD2*ROQG)
       EGR=1.
      RN16=CPI2*EGR*TNG*TNW*ROQR
      RN17=2.*CPI*TNG
       RN17A=.31*GA5GH*(CD2*ROQG)**.25
       RN17B=CW-CI
       RN17C=CW
       APRI=.66
       BPRI=1.E-4
      RN18=20.*CPI2*BPRI*TNW*ROQR
       RN18A=APRI
      RN19=2.*CPI*TNG/ALF
       RN19A=.31*GA5GH*(CD2*ROQG)**.25
       RN19B=CW/ALF
      RN20=2.*CPI*TNG
       RN20A=ALS*ALS/RW
       RN20B=.31*GA5GH*(CD2*ROQG)**.25
       BND3=2.E-3
      RN21=1.E3*1.569E-12/0.15
       ERW=1.
      RN22=.25*CPI*ERW*AC1*TNW*GA3B
      RN23=2.*CPI*TNW
       RN23A=.31*GA5BH*SQRT(AC1)
       RN23B=ALV*ALV/RW
       CN0=1.E-8
      RN25=CN0
       RN30A=ALV*ALS*AMW/(TCA*ARS)
       RN30B=ALV/TCA
       RN30C=ARS/(DWV*AMW)
      RN31=1.E-17
       BETA=-.6
      RN32=4.*51.545E-4
c     ****************************
      do 20 k=1,nz
       srr(k)=1./sqrt(rho(k))
   20  qrr(k)=sqrt(srr(k))
      RETURN
      END
C-----------------------------------------------------------------------------
      subroutine base (IRS)
      PARAMETER (NX=514,NZ=43,ITT=244,NX10=10*NX)
      common/msui/ isounding,isui,ifine,idelz
      COMMON/OPTION/ LIPPS,IMLIFTING,ISNGTAO,IJKADV,IWBAR,iuvbar
      common/iceopt/ ice913,ilif
      common/cpsbm/ icps_bm,iexplicit
      common/o4/ a4k(nz),a4h(nz),d58x,d16x,d48x,d516x,d32x,d96x
      common/bxz/ imax,iles,il2,kmax,kles,kl2
      common/bx/ dx,dx2,rdx,rd2x,rd4x,rdx2,r2dx2,r4dx2
      common/bz/ dz,dz2,rdz,rd2z,rd4z,rdz2,r2dz2,r4dz2
      common/bb/ dt,d2t,ril2,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,psfc
      common/b1u/ u(nx,nz)
      common/b1v/ v(nx,nz)
      common/b1w/ w(nx,nz)
      common/b2u/ uu1(nx,nz)
      common/b2v/ vv1(nx,nz)
      common/b2w/ ww1(nx,nz)
      common/dumuw/ umd(nx,nz),vmd(nx,nz),wmd(nx,nz)
      common/b4/ tbsk(nz),bskt(nz),bskt2(nz),bskm(nz),bskm4(nz),
     $   bsit(nz),bsit2(nz),bsim(nz),bsim4(nz)
      common/b5/ tb(nz),qb(nz),rho1(nz),rho(nz),ba(nz),bb(nz),ta(nz),
     $   qa(nz),ta1(nz),qa1(nz),coef(nz),c1(nz),c2(nz),c3(nz),am(nz),
     $   am1(nz),ub(nz),vb(nz),wb(nz),ub1(nz),vb1(nz),rrho(nz),rrho1(nz)
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     $   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
      common/bb6/ tls1(nz),tls2(nz),qls1(nz),qls2(nz),tls3(nz),tls4(nz),
     $  qls3(nz),qls4(nz),sft(nz),sfq(nz),wbt(nz),wb_6h(nz,itt),
     $  ub_6h(nz,itt),ubt(nz),q1_6h(nz,itt),q1t(nz),q2_6h(nz,itt),
     *  q2t(nz),vb_6h(nz,itt),vbt(nz)

      COMMON/UPPER1/ T_ADJUST(NZ,ITT),Q_ADJUST(NZ,ITT),PRESS_T(NZ),
     1     TEMP_T(NZ),VAP_T(NZ)
      common/bls3/ factor_nuding
      common/q1q2t/ aq1t(nz),aq2t(nz),aq1zt(nz),aq2zt(nz)
      common/q1q2z/ Q1Z_6H(NZ,ITT),Q1ZT(NZ),Q2Z_6H(NZ,ITT),Q2ZT(NZ)
      common/ba/ y1(nx),y2(nx),y3(nx),qvs(nx),qr(nx),hu1(nx),y4(nx10)
      common/damp/ rfa(nz),rfa1(nz),cntc(nz),cgwd
      common/bls/ wbx(nx),ths(nx,itt),qs(nx,itt),ts(nx,itt),pss(nx,itt)
      common/dinrad/p00(nz),dz0(nz),tairsfc(nx),qairsfc(nx),thairsf(nx)
     *,       pairsfc(nx)
      common/dinrad1/ dz1half
      common/bpbl/ UHT(NZ),WHT(NZ),TGBAT0
      real z1(nz),z2(nz),hgt1(nz),hgt2(nz)
      common/zlevel/ cc1,cc2,z1,z2
      common/bnumer/ kref
      real dzy(nz),dzyp(nz),sst(itt),ssp(itt)
      save
c      DATA DZY/0.,120.,150.,180.,210.,240.,270.,300.,350.,400.,450.,
c     1   500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,
c     2   1000.,1000.,1000.,1000.,1000.,1100.,1100.,1100.,1100.,1100.,
c     3   1200.,1200./
c      DATA DZYP/0.,120.,150.,180.,210.,240.,270.,300.,350.,400.,450.,
c     1   500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,
c     2   1000.,1000.,1000.,1000.,1000.,1100.,1100.,1100.,1100.,1100.,
c     3   1200.,1200./

      ILES=NX-2
      KLES=NZ-2
      IL2=ILES
      KL2=KLES
      ILES=ILES+1
      KLES=KLES+1
      IMAX=ILES+1
      KMAX=KLES+1
      RIL2=1./FLOAT(IL2)

C Grid setup for GCSS WG4 model intercomparison workshop

      dzy(1)=0.
      dzyp(1)=0.
      dzy(2)=0.
      dzyp(2)=0.
      do k=3,nz
         if(k.le.22)then
            dzk=real(k-2)*300.
            dzkp=real(k-2.5)*300.
            dzy(k)=(.2982456142+1.169590643E-4*dzk)*dzk
            dzyp(k)=(.2982456142+1.169590643E-4*dzkp)*dzkp
         else
            dzy(k)=real(k-22)*500.+6000.
            dzyp(k)=real(k-22.5)*500.+6000.
         endif
      enddo
      do k=2,nz-1
         dzyp(k)=dzyp(k+1)
      enddo
      dzyp(nz)=dzyp(nz-1)+500.

C End GCSS WG4 grid setup

      if(idelz.eq.1)then 
         DZ=DZY(3)*100.
      else

C  IF USING 34 PTS, SET DZ = 27.e2
C  IF USING 40 PTS, SET DZ = 25.e2
C  IF USING 43 PTS, SET DZ = 20.e2

         DZ=20.e2      !cm

      endif
cccshie
      DX=1000.e2    !cm for nx= 514, 5/18-5/26/98
c     DX= 500.e2    !cm for nx= 1026, 6/2-6/11/98
c     DX=2000.e2    !cm for nx= 258, 5/18-5/20/98
c
      IF(ICPS_BM .EQ. 1) DX=30000.E2
c
      PSFC=1007.e3
      TSFC=27.2     !C

C IF USING 34 PTS, SET CC1=4, CC2=1/40.
C IF USING 40 PTS, SET CC1=6, CC2=1/50.
C IF USING 43 PTS, SET CC1=3.4, CC2=.03

      CC1=3.4
      CC2=.03

      RD1=1.E-3
      RD2=2.2
      BOUND=1.E-3
      AL=2.5E10
      CP=1.004E7
      RA=2.87E6
      F5=17.26939*237.3*AL/CP
      SCALE=-1.*1000./3600.

      Z1(1) = 0.
      Z2(1) = 0.

      if (idelz.eq.1) then
         do k=1,kmax
            Z2(k)=dzy(k)*100.
            Z1(k)=dzyp(k)*100.
         enddo
      else
         DO K=2,KMAX
            if (ifine .eq. 0) then
               Z1(k)=(cc1+cc2*(k-1.5)*dz*.01)*(k-1.5)*dz
               Z2(k)=(cc1+cc2*(k-2.)*dz*.01)*(k-2.)*dz
            else
               y1(k)=(cc1+cc2*(k-1.5)*dz*.01)*(k-1.5)*dz
               y2(k)=(cc1+cc2*(k-2.)*dz*.01)*(k-2.)*dz
            endif

         ENDDO
         if (ifine .eq. 1) then
            Z1(2)=40.E2
            Z2(2)=0.0
            Z2(3)=80.E2
            DO K=2,KLES
               Z1(K+1)=Y1(K)+80.E2
               Z2(K+1)=Y2(K)+80.E2
               print*,k,z2(k)/100.,(z2(k)-z2(k-1))/100.
            ENDDO
         endif
      endif

      DO K=2,KLES
         AM(K)=DZ/(Z2(K+1)-Z2(K))
         AM1(K)=DZ/(Z1(K)-Z1(K-1))
      ENDDO
      AM1(KMAX)=DZ/(Z1(KMAX)-Z1(KLES))
      AM1(1)=AM1(2)
      AM(1)=AM(2)
      AM(KMAX)=AM(KLES)
      do k=1,kmax
        hgt1(k)=z1(k)
        hgt2(k)=z2(k)
        uht(k)=z1(k)
        wht(k)=z2(k)
      enddo
      dz1half=z1(2)/100.

      write(6,*)' km      am        am1        z1          z2       
     1dz2       dratio'
      do 14 k=kmax,2,-1
      deltz2=z2(k)-z2(k-1)
      if(k.ge.4)then
         deltz3=z2(k-1)-z2(k-2)
         dratio=deltz2/deltz3
      else
         dratio=0.
      endif
   14 write(6,2035) k,am(k),am1(k),z1(k)/100.,z2(k)/100.,
     1              deltz2/100.,dratio
      write(6,*)
      write(6,*)

      DX2=DX*DX
      DZ2=DZ*DZ
      RDX=1./DX
      RDZ=1./DZ
      RD2X=.5*RDX
      RD2Z=.5*RDZ
      RD4X=.25*RDX
      RD4Z=.25*RDZ
      RDX2=RDX*RDX
      RDZ2=RDZ*RDZ
      R2DX2=.5*RDX2
      R2DZ2=.5*RDZ2
      R4DX2=.25*RDX2
      R4DZ2=.25*RDZ2
      D16X=1./16.*RDX
      D48X=1./48.*RDX
      D58X=5./8.*RDX
      D32X=1./32.*RDX
      D96X=1./96.*RDX
      D516X=5./16.*RDX

C     ******   RAYLEIGH RELAXATION TERMS 
      RRTP0=10.
      RRTP=15.5  ! tao 5/1/01
c     RRTP=16.
      RRTPD=RRTP-RRTP0
      CGWD=1./(5.*24.*3600.)
      CG=1./(24.*3600.)
      CRF=1./900.
      DO 15 K=1,KMAX
        A1=1.E-5*Z1(K)
        A2=1.E-5*Z2(K)
        UHT(K)=.01*Z1(K)
        WHT(K)=.01*Z2(K)
        RFA(K)=0.0
        RFA1(K)=0.0
        CNTC(K)=0.0
      IF(A1.GE.RRTP) CNTC(K)=CG
      IF(A1.LT.RRTP .AND. A1.GT.RRTP0) CNTC(K)=CG*(A1-RRTP0)/RRTPD
      IF(A1.GE.RRTP) RFA(K)=CRF*(A1-RRTP)
   15 IF(A2.GE.RRTP) RFA1(K)=CRF*(A2-RRTP)
      CKH=2.
      CKHB=2.
c
      A2K=0.2E6
c
      AAK=0.00175/DT
c
      if (ijkadv .eq. 1) then
         A2K=0.15E6 ! tao 5/1/01
c        A2K=0.2E6
c         aak=0.00175/dt	!tao-new
         aak=0.001250/dt  ! tao 5/1/01
      endif

      IF(ICPS_BM .EQ. 1) a2k=.1e6
      IF(ICPS_BM .EQ. 1) AAK=0.00125/DT
c
      kref=2
      zkref=250.    !m
      do k=2,kmax
        if(((z2(k)-z2(k-1))/100.).le.zkref) kref=k
      enddo
      print*,'kref is',kref
      A2I=AAK*DX*DX
      DO 20 K=1,KMAX
         Y3(K)=DZ*DZ/(AM(K)*AM(K))
         COEF(K)=CK*CK*Y3(K)*.5
         C1(K)=3.*980.*COEF(K)
         C2(K)=.5*CE/(CK*Y3(K))
         C3(K)=2./(3.*CK*CK*Y3(K))
ctao  (11-13)
        A4H(K)=AAK*(AM(kref)/AM(K))
c        A4H(K)=AAK
          IF(ICPS_BM .EQ. 1) A4H(K)=AAK*(AM(kref)/AM(K))
          A4K(K)=A4H(K)/CKHB
c
        TBSK(K)=2.0e6*(AM(kref)/AM(K)) ! tao 5/1/01
c       TBSK(K)=2.5e6*(AM(kref)/AM(K))
         IF(ICPS_BM .EQ. 1) TBSK(K)=2.5e6*(AM(kref)/AM(K))/10.
c
        BSKT(K)=A2K*(AM(kref)/AM(K))
          BSKT2(K)=2.*BSKT(K)
         BSKM(K)=BSKT(K)/CKH
           BSKM4(K)=4.*BSKM(K)
c
         BSIT(K)=A2I*(AM(kref)/AM(K))
           BSIT2(K)=2.*BSIT(K)
          BSIM(K)=BSIT(K)/CKH
            BSIM4(K)=4.*BSIM(K)
   20 CONTINUE
      if(isounding.eq.1)then 
        READ(5,202) (FD(K),TB(K),K=2,KMAX)
        READ(5,202) (FE(K),QB(K),K=2,KMAX)
        READ(5,250) (UB(K),K=2,KMAX)
        READ(5,250) (VB(K),K=2,KMAX)
        READ(5,202) (WB(K),K=2,KMAX)
        IF (IWBAR .EQ. 1) THEN
          DO ICY=1,itt
            READ(5,202) (WB_6H(K,ICY),K=2,KMAX)
          ENDDO
          IF (ISUI .EQ. 0) THEN
            DO ICY=1,itt
              READ(5,202) (Q1Z_6H(K,ICY),K=2,KMAX)
              READ(5,202) (Q2Z_6H(K,ICY),K=2,KMAX)
            ENDDO
          ENDIF
          DO K=2,KMAX
            QR(K)=WB(K)
            DO ICY=1,itt
              QR(K)=QR(K)+WB_6H(K,ICY)
            ENDDO
            QR(K)=QR(K)/(itt+1)
            WB(K)=QR(K)
          ENDDO
        ELSE
          DO ICY=1,itt
            READ(5,250) (UB_6H(K,ICY),K=2,KMAX)
            READ(5,250) (VB_6H(K,ICY),K=2,KMAX)
            READ(5,202) (WB_6H(K,ICY),K=2,KMAX)
          ENDDO
          do k=23,kmax
            ub(k)=ub(22)
            do icy=1,itt
              ub_6h(k,icy)=ub_6h(22,icy)
            enddo
          enddo
          DO ICY=1,itt
            READ(5,202) (Q1_6H(K,ICY),K=2,KMAX)
            READ(5,202) (Q1Z_6H(K,ICY),K=2,KMAX)
            READ(5,202) (Q2_6H(K,ICY),K=2,KMAX)
            READ(5,202) (Q2Z_6H(K,ICY),K=2,KMAX)
          ENDDO
          DO ICY=1,itt
            Q1Z_6H(1,ICY)=Q1Z_6H(2,ICY)
            Q2Z_6H(1,ICY)=Q2Z_6H(2,ICY)
          ENDDO
        ENDIF
      else
c
        call ainter(sst,ssp,tb,fd,qb,fe,ub,vb,wb,z1,z2,toptq,topuv)
c
        PSFC=SSP(1)*1000.
        do k=1,kmax
          if((z2(k)/100.).le.toptq) ktoptq=k
          if((z2(k)/100.).le.topuv) ktopuv=k
c          do i=1,itt
c            q1z_6h(k,i)=0.
c            q2z_6h(k,i)=0.
c          enddo
        enddo
       print*,'ktoptq, m, tophtq ',ktoptq,.01*z2(ktoptq),toptq
       print*,'ktopuv, m, tophuv ',ktopuv,.01*z2(ktopuv),topuv
      endif
c      rate=5.
c      rate=max(rate,tb(kmax)-tb(ktoptq))
c      rate=min(rate,20.)/(real(kmax-ktoptq)+1.e-10)
c      do k=ktoptq+1,kmax
c        tb(k)=rate+tb(k-1)
c        fd(k)=rate+fd(k-1)
c        qb(k)=0.
c        fe(k)=0.
c        wb(k)=0.
c      enddo
c
cc
ccc

cc
c
      IF (ILIF .EQ. 1) THEN
c       dili is used to read v,tvad,qvad
        r91=1.0
        r92=0.1
        z91=z1(18)
        z92=z1(26)
        a99=(r91-r92)/(z91-z92)
        b99=(r91*z92-r92*z91)/(z92-z91)
        do k=1,9
          z99=z1(17+k)
          r99=a99*z99+b99
         if(k.eq.1) ry1=r99
         if(k.eq.2) ry2=r99
         if(k.eq.3) ry3=r99
         if(k.eq.4) ry4=r99
         if(k.eq.5) ry5=r99
         if(k.eq.6) ry6=r99
         if(k.eq.7) ry7=r99
         if(k.eq.8) ry8=r99
         if(k.eq.9) ry9=r99
       enddo
       write(6,*) 'test',ry1,ry2,ry3,ry4,ry5,ry6,ry7,ry8,ry9
        wb(18)=ry1*wb(18)
	wb(19)=ry2*wb(19)
	wb(20)=ry3*wb(20)
	wb(21)=ry4*wb(21)
	wb(22)=ry5*wb(22)
	wb(23)=ry6*wb(23)
	wb(24)=ry7*wb(24)
	wb(25)=ry8*wb(25)
	wb(26)=ry9*wb(26)
	wb(27)=ry9*wb(27)
	wb(28)=ry9*wb(28)
	wb(29)=ry9*wb(29)
	wb(30)=ry9*wb(30)
	wb(31)=ry9*wb(31)
	wb(32)=ry9*wb(32)
	wb(33)=ry9*wb(33)
        do icy=1,ITT
           wb_6h(18,icy)=ry1*wb_6h(18,icy)
           wb_6h(19,icy)=ry2*wb_6h(19,icy)
	   wb_6h(20,icy)=ry3*wb_6h(20,icy)
	   wb_6h(21,icy)=ry4*wb_6h(21,icy)
 	   wb_6h(22,icy)=ry5*wb_6h(22,icy)
	   wb_6h(23,icy)=ry6*wb_6h(23,icy)
	   wb_6h(24,icy)=ry7*wb_6h(24,icy)
	   wb_6h(25,icy)=ry8*wb_6h(25,icy)
	   wb_6h(26,icy)=ry9*wb_6h(26,icy)
	   wb_6h(27,icy)=ry9*wb_6h(27,icy)
	   wb_6h(28,icy)=ry9*wb_6h(28,icy)
	   wb_6h(29,icy)=ry9*wb_6h(29,icy)
	   wb_6h(30,icy)=ry9*wb_6h(30,icy)
	   wb_6h(31,icy)=ry9*wb_6h(31,icy)
	   wb_6h(32,icy)=ry9*wb_6h(32,icy)
	   wb_6h(33,icy)=ry9*wb_6h(33,icy)
        enddo
        do k=23,kmax
          ub(k)=ub(22)
        do icy=1,ITT
          ub_6h(k,icy)=ub_6h(22,icy)
        enddo
        enddo
        do icy=1,ITT
	   q1_6h(18,icy)=ry1*q1_6h(18,icy)
	   q1_6h(19,icy)=ry2*q1_6h(19,icy)
	   q1_6h(20,icy)=ry3*q1_6h(20,icy)
	   q1_6h(21,icy)=ry4*q1_6h(21,icy)
	   q1_6h(22,icy)=ry5*q1_6h(22,icy)
	   q1_6h(23,icy)=ry6*q1_6h(23,icy)
	   q1_6h(24,icy)=ry7*q1_6h(24,icy)
	   q1_6h(25,icy)=ry8*q1_6h(25,icy)
	   q1_6h(26,icy)=ry9*q1_6h(26,icy)
	   q1_6h(27,icy)=ry9*q1_6h(27,icy)
	   q1_6h(28,icy)=ry9*q1_6h(28,icy)
	   q1_6h(29,icy)=ry9*q1_6h(29,icy)
	   q1_6h(30,icy)=ry9*q1_6h(30,icy)
	   q1_6h(31,icy)=ry9*q1_6h(31,icy)
	   q1_6h(32,icy)=ry9*q1_6h(32,icy)
	   q1_6h(33,icy)=ry9*q1_6h(33,icy)
cxx    do k=2,kmax
cxx       q1_6h(k,icy)=q1_6h(k,icy)+273.16
cxx       q2_6h(k,icy)=q2_6h(k,icy)*0.001
cxx    enddo
       enddo
       wb(1)=0.
       wb(2)=0.
       wb(kmax)=0.
       do icy=1,ITT
          ub_6h(   1,icy)=ub_6h(   2,icy)
          vb_6h(   1,icy)=vb_6h(   2,icy)
          q1_6h(   1,icy)=q1_6h(   2,icy)
          q2_6h(   1,icy)=q2_6h(   2,icy)
          wb_6h(   1,icy)=0.
          wb_6h(   2,icy)=0.
          wb_6h(kmax,icy)=0.
       end do
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,210)
      DO 32 K1=2,KMAX
        K=KMAX+2-K1
   32 WRITE(6,2031) K,TB(K),QB(K),FD(K),FE(K),QR(K),UB(K),VB(K),WB(K)
      WRITE(6,204)
      DO 34 K1=2,KMAX
        K=KMAX+2-K1
   34 WRITE(6,2441) K,Z1(K)/100.,Z2(K)/100.,AM(K),AM1(K),COEF(K)
     *,              C1(K),C2(K),C3(K),TBSK(K),BSKM(K),RFA(K)
      DO 40 K=2,KMAX
        FD(K)=FD(K)+273.16
        FE(K)=FE(K)*1.E-3
        TB(K)=TB(K)+273.16
   40   QB(K)=QB(K)*1.E-3
      A1=980.*DZ/RA
      P0(2)=LOG(PSFC)
      DO 45 K=3,KMAX
   45   P0(K)=P0(K-1)-A1/(TB(K-1)*(1.+.61*QB(K-1))*AM(K-1))

      RHO1(KMAX)=EXP(P0(KMAX))/(RA*FD(KMAX)*(1.+.61*FE(KMAX)))
      DO 46 K=2,KLES
   46   DZ0(K)=Z2(K+1)-Z2(K)
      DZ0(KMAX)=DZ0(KLES)
      DO 47 K=2,KMAX
   47   P00(K)=EXP(P0(K))
      DO 50 K=2,KLES
        RHO1(K)=EXP(P0(K))/(RA*FD(K)*(1.+.61*FE(K)))
        P0(K)=EXP((P0(K)+P0(K+1))*.5)
        PI(K)=(P0(K)/1000.E3)**.286
        F0(K)=AL/(PI(K)*CP)
        Y1(K)=TB(K)
        TB(K)=TB(K)/PI(K)
        QVS(K)=3.799E3/P0(K)*EXP(17.26939*(Y1(K)-273.16)/(Y1(K)-35.86))
        Y3(K)=3.799E3/P0(K)*EXP(21.87456-5807.695/(Y1(K)-7.66))
c        if(isounding.eq.1)then
c          if (k.eq.26) qb(k)=0.3*y3(k)
c          if (k.eq.27) qb(k)=0.225*y3(k)
c          if (k.eq.28) qb(k)=0.15*y3(k)
c          if (k.eq.29) qb(k)=0.015*y3(k)
c          if (k.eq.30) qb(k)=0.005*y3(k)
c          if (k.eq.31) qb(k)=0.001*y3(k)
c          if (k.ge.32) qb(k)=0.0   
c        endif
        Y4(K)=QB(K)/Y3(K)
        HU1(K)=QB(K)/QVS(K)
c tao add 5/01/01, try to resolve over-cooling problem
        if (k .ge. 23 .and. y4(k) .ge. 0.50) qb(k)=0.5*y3(k)
          if (k.eq.24) qb(k)=0.45*y3(k)
          if (k.eq.25) qb(k)=0.40*y3(k)
          if (k.eq.26) qb(k)=0.35*y3(k)
          if (k.eq.27) qb(k)=0.30*y3(k)
          if (k.eq.28) qb(k)=0.25*y3(k)
          if (k.eq.29) qb(k)=0.15*y3(k)
          if (k.eq.30) qb(k)=0.10*y3(k)
          if (k.eq.31) qb(k)=0.075*y3(k)
          if (k.ge.32) qb(k)=0.025*y3(k)  
          if (k.ge.33) qb(k)=0.0   
        Y4(K)=QB(K)/Y3(K)
        HU1(K)=QB(K)/QVS(K)

   50 RHO(K)=P0(K)/(RA*Y1(K)*(1.+.61*QB(K)))
      DO K=3,KLES-1
        IF(TB(K).LT.TB(K-1)) TB(K)=.5*(TB(K-1)+TB(K+1))
        IF(TB(K).LT.TB(K-1)) TB(K)=TB(K-1)-TB(K)+TB(K-1)
      ENDDO
      IF(TB(KLES).LT.TB(KLES-1))TB(KLES)=TB(KLES-1)*2.-TB(KLES-2)
        rqm=0.
        do k=2,kles
        rqm=rqm+(.5*(rho1(k)+rho(k))*.5*(qb(k)+fe(k)))*(hgt1(k)-hgt2(k))
        rqm=rqm+(.5*(rho1(k+1)+rho(k))*.5*(qb(k)+fe(k+1)))
     *     *(hgt2(k+1)-hgt1(k))
        enddo 
        print*,'MODEL QV INTEGRATION=',rqm
        print*

        if (ilif .eq. 0) then

         zpresst=75.*1000.
         zpressu=75.*1000.
        do k=2,kles
           if(p0(k) .ge. zpresst) toptq=k
           if(p0(k) .ge. zpressu) topuv=k
        enddo
        print*,'top at or just below 75mb =',toptq
        print*,'top at or just below 75mb =',topuv
          do k=ktopuv+1,kmax
            ub(k)=ub(ktopuv)
            vb(k)=vb(ktopuv)
          enddo
         do k=ktoptq+1,kmax
           wb(k)=0.
         enddo

cccshie  3/30/01 by dan: ADJUST THE VALUES NEAR THE TOP

         do k=2,kles
           p0mb=p0(k)*.001 
c          if (p0mb.lt.350.) then
           if (p0mb.lt.450.) then  ! cccshie by tao 5/3/01
              k350=k
              goto 4321
           endif
         enddo
4321     continue

         do k=2,kles
           p0mb=p0(k)*.001 
c          if (p0mb.lt.100.) then
           if (p0mb.lt.150.) then ! cccshie by tao 5/4/01
              k100=k
              goto 4323
           endif
         enddo
4323     continue

cccshie by tao 5/4/01
         do k=2,kles
           p0mb=p0(k)*.001
           if (p0mb.lt.250.) then
              k250=k
              goto 4325
           endif
         enddo
4325     continue

         do icy=1,itt
            do k=ktopuv+1,kmax
              ub_6h(k,icy)=ub_6h(ktopuv,icy)
              vb_6h(k,icy)=vb_6h(ktopuv,icy)
            enddo
            do k=ktoptq+1,kmax
              wb(k)=0.
              wb_6h(k,icy)=0.
c              q1_6h(k,icy)=0.   ! dan 3/30/01
c              q2_6h(k,icy)=0.   ! dan 3/30/01
c              q1z_6h(k,icy)=0.  ! dan 3/30/01
c              q2z_6h(k,icy)=0.  ! dan 3/30/01
            enddo
            wb_6h(2,icy)=0.
            wb_6h(1,icy)=0.
         
cccshie  3/30/01 by dan: ADJUST THE VALUES NEAR THE TOP
         do k=k350,k100
c           divfac=(1.-real(k-k350)/real(k100-k350))
            divfac=(1.-real(k-k350)/real(k100-k350))**2. ! cccshie by tao 5/3/01
            q1_6h(k,icy)=q1_6h(k350,icy)*divfac
c           q2_6h(k,icy)=q2_6h(k350,icy)*divfac    ! cccshie by tao 5/4/01
            q1z_6h(k,icy)=q1z_6h(k350,icy)*divfac
c           q2z_6h(k,icy)=q2z_6h(k350,icy)*divfac  ! cccshie by tao 5/4/01
         enddo
cccshie 5/4/01 shie modified following tao's idea: to drier top layers sounding
         do k=k350,k250
           divfacv=(1.-real(k-k350)/real(k250-k350))**2.
            q2_6h(k,icy)=q2_6h(k350,icy)*divfacv
            q2z_6h(k,icy)=q2z_6h(k350,icy)*divfacv
         enddo
         do k=k250,k100
            q2_6h(k,icy)=0.
            q2z_6h(k,icy)=0.
         enddo

         do k=k100+1,kmax
            q1_6h(k,icy)=0.
            q2_6h(k,icy)=0.
            q1z_6h(k,icy)=0.
            q2z_6h(k,icy)=0.
         enddo

       enddo
      endif

cc
cc
      IF (LIPPS.EQ.0)THEN
        DO K=3,KLES
          A1=RHO1(K)*.5
          FD(K)=(TB(K)-TB(K-1))*A1
          FE(K)=(QB(K)-QB(K-1))*A1
          A2=AM1(K)*(TB(K)+TB(K-1))*A1
          BA(K)=AM(K)*A2/(RHO(K)*TB(K))
          BB(K-1)=AM(K-1)*A2/(RHO(K-1)*TB(K-1))
        ENDDO
      ELSE
        DO K=3,KLES
          A1=RHO1(K)*.5
          FD(K)=(TB(K)-TB(K-1))*A1
          FE(K)=(QB(K)-QB(K-1))*A1
          A2=AM1(K)*RHO1(K)
          BA(K)=AM(K)*A2/RHO(K)
          BB(K-1)=AM(K-1)*A2/RHO(K-1)
        ENDDO
      ENDIF
      BA(2)=AM(2)*AM1(2)*RHO1(2)/RHO(2)
      BB(KLES)=AM(KLES)*AM1(KMAX)*RHO1(KMAX)/RHO(KLES)
c      if(iuvbar.eq.0)then 
c        A2=0.
c        A3=0.
c        A4=0.
c        DO 6666 K=2,KLES
c          A1=RHO(K)*DZ/AM(K)
c          A2=A2+A1*UB(K)
c          A3=A3+A1*VB(K)
c          A4=A4+A1
c 6666   CONTINUE
c        DO 6669 K=2,KLES
c 6669     VB(K)=A3/A4
c      endif
      do k=itt,1,-1
        if(isounding.eq.1)then
          TSL=(TSFC+273.16)/pi(2)
          QSL=3.799052E3/p0(2)*EXP(17.26939*TSFC/(TSFC+237.3))*.975
        else 
          TSFC=SST(k)
          PSFC=SSP(k)*1000.
          pisfc=(psfc/1000.e3)**.286
          TSL=(TSFC+273.16)/pisfc
          QSL=3.799052E3/psfc*EXP(17.26939*TSFC/(TSFC+237.3))
        endif
          print*,'psfc',psfc,psfc*.001
        write(6,12121) k,tsl,qsl*1000.,psfc*.001,tsfc+273.16
        do i=1,imax
          THS(i,k)=TSL
          QS(i,k)=QSL
          TS(i,k)=tsfc+273.16
          PSS(i,k)=psfc
        enddo
      enddo
12121 format(2x,3hit=,i4,4x,5htsfc=,f12.5,4x,5hqsfc=,f12.5,4x,
     1          5hpsfc=,f12.5,4x,7htsfc c=,f12.5)

      TB(1)=TB(2)
      QB(1)=QB(2)
      UB(1)=UB(2)
      VB(1)=VB(2)
      WB(1)=0.
      wb(kmax)=0.
      PI(1)=PI(2)
      HU1(1)=HU1(2)
      FD(2)=0.
      FE(2)=0.
      FD(1)=FD(2)
      FE(1)=FE(2)
      F0(1)=F0(2)
      RHO1(1)=RHO1(2)
      RHO(1)=RHO(2)
      P0(1)=P0(2)
      TB(KMAX)=TB(KLES)*2.-TB(KL2)
      QB(KMAX)=QB(KLES)*2.-QB(KL2)
      UB(KMAX)=UB(KLES)
      VB(KMAX)=VB(KLES)
      PI(KMAX)=PI(KLES)
      HU1(KMAX)=0.
      F0(KMAX)=2.*F0(KLES)-F0(KL2)
      FD(KMAX)=FD(KLES)
      FE(KMAX)=FE(KLES)
      RHO(KMAX)=RHO(KLES)
      P0(KMAX)=P0(KLES)
      DO 70 K=1,KMAX
        wb(k)=scale*wb(k)/(rho1(k)*980.)
        if(k.le.2) wb(k)=0.
        if(k.eq.kmax) wb(k)=0.
           if (p0(k) .le. 75.e3) then
             if (wb(k) .lt. 0.0) wb(k)=0.
           endif
         do icy=1,itt
          if (k.le.2) wb_6h(k,icy)=0.
          if (k.eq.kmax) wb_6h(k,icy)=0.
           wb_6h(k,icy)=scale*wb_6h(k,icy)/(rho1(k)*980.)
             if (p0(k) .le. 75.e3) then
               if (wb_6h(k,icy) .lt. 0.0) wb_6h(k,icy)=0.
             endif
         enddo
        RRHO(K)=1./RHO(K)
        RRHO1(K)=1./RHO1(K)
   70 CONTINUE
      IF(IRS.EQ.1) GO TO 1111
      DO 150 K=1,KMAX
        TA(K)=TB(K)
        QA(K)=QB(K)
        TA1(K)=TB(K)
        QA1(K)=QB(K)
        UB1(K)=UB(K)
        VB1(K)=VB(K)
      DO 150 I=1,IMAX
        UU1(I,K)=UB(K)
        VV1(I,K)=VB(K)
        ww1(i,k)=0.
        U(I,K)=UB(K)
        V(I,K)=VB(K)
        w(i,k)=0.
        umd(I,K)=UB(K)
        vmd(i,k)=vb(k)
        wmd(I,K)=0.
  150 CONTINUE
 1111 WRITE(6,200)
      DO 33 K1=2,KLES
      K=KLES+2-K1
      KM=K-1
   33 WRITE(6,2032)KM,P0(K)/1000.,PI(K),TB(K),QB(K)*1000.,QVS(K)*1000.
     *,     HU1(K),BA(K)*1000.,BB(K)*1000.,Y3(K)*1000.,Y4(K)
      WRITE(6,201)
      DO 35 K1=2,KLES
      K=KLES+2-K1
      KM=K-1
   35 WRITE(6,2033)KM,RHO(K)*1000.,RHO1(K)*1000.,F0(K),FD(K),FE(K),
     1             UB(K)/100.,VB(K)/100.,WB(K),BSIM(K)
      do icy=1,itt
        print*
        print*,'itt =',icy
        print*
        write(6,241)
        do 36 k1=2,kles
          k=kles+2-k1
          km=k-1
   36     write(6,203) km,P0(K)/1000.,ub_6h(k,icy),vb_6h(k,icy),
     1                   wb_6h(k,icy),Q1_6H(k,icy),Q2_6H(k,icy)*1000.,
     2                   t_adjust(k,icy),q_adjust(k,icy)*1000.
      enddo
      RETURN
  200 FORMAT(//,' LEVEL   PRES    PI      TB       QB     QSW     HUM      
     *   BA      BB     QSI      HUI')
  201 FORMAT(//,'LEVEL  RHO   RHO1    F0        FD       FE       UB   
     1     VB    WB      BSIM')
  202 FORMAT(16F5.2)
  203 FORMAT(I4,10F12.5,/)
 2035 FORMAT(I4,2F10.5,4F12.3)
 2031 FORMAT(1X,I4,2(F8.2,F7.2),F8.2,2F10.2,F8.2)
 2032 FORMAT(I4,F9.2,F7.3,F9.2,7F8.2)
 2033 FORMAT(I4,2F7.3,F9.2,E11.3,E12.3,2F7.2,F7.3,E11.3)
  204 FORMAT(//,'LEVEL  Z      Z1     AM    AM1    COEF       C1       
     *C2       C3        TBSK       BSKM')
 2441 FORMAT(I3,2F8.1,2F6.3,5E10.3,2E9.2,E10.3)
  210 FORMAT(//,'  LEVEL  TB      QB     FD     FE1      QR       UB   
     1     VB        WB')
  241 format(//,'level          p           u           v           w 
     1       q1          q2        tdt       qdt')
  250 format(16f5.0)
      END
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

C --------------------------------------------------------------------------

      FUNCTION RAN3 (idum)
      implicit none
      integer mbig,mseed,mz,idum
      real    ran3,fac
      PARAMETER  (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
      integer i,k,iff,mj,mk,ii,inext,inextp,ma(55)
      SAVE  inext,inextp,ma,iff
      DATA  iff /0/
      IF( idum.LT.0 .OR. iff .EQ. 0 ) THEN
        iff= 1
        mj= mseed - ABS( idum )
c       mj=MOD( mj,mbig )
        mj= mj-mbig*(mj/mbig)
        ma(55)= mj
        mk= 1
        DO 10 i=1,54
c          ii= MOD( 21*i,55 )
          ii= (21*i)-55*((21*i)/55)
          ma(ii)= mk
          mk= mj-mk
          IF( mk.LT.mz ) mk= mk+mbig
          mj= ma(ii)
 10     CONTINUE
        DO 20 k=1,4
          DO 15 i=1,55
c            ma(i)= ma(i) - ma(1+MOD( i+30,55 ))
            ma(i)= ma(i) - ma(1+(i+30)-55*((i+30)/55))
            IF( ma(i).LT.mz ) ma(i)= ma(i)+mbig
 15       CONTINUE
 20     CONTINUE
        inext= 0
        inextp= 31
        idum= 1
      ENDIF
      inext= inext+1 
        IF( inext .EQ. 56 ) inext= 1
      inextp= inextp+1
        IF( inextp .EQ. 56 ) inextp= 1
      mj= ma(inext) - ma(inextp)
        IF( mj .LT. mz ) mj= mj+mbig
      ma(inext)= mj
      ran3= REAL( mj )*fac
      RETURN

      END
C  FILE TO READ IN OBS. DATA FOR GATE.  FILES ARE ORGINALLY IN BINARY FORMAT

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION CVMGP(X1,X2,X3)
       IF(X3.GE.0.)  THEN
        CVMGP=X1
       ELSE
        CVMGP=x2
       ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine budget_m_rh (ijkadv,ndt_stat,d2t)
      parameter (NX=514,NZ=43)
c     modified for TOGA COARE and GATE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/mbudget/ acoc(nz),acoe(nz),agaut(nz),asdep(nz),agsub(nz),
     1   apern(nz),aqls(nz),altrans(nz),actrans(nz),total_mq(nz),
     2   agwet(nz),afgr(nz),avelw(nz),avesw(nz),amelt(nz),aother(nz),
     3   total_mt(nz)
      common/tqave/ tavet(nz),qavet(nz)

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

      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNQH(NZ,4,7),
     1 SNQV(NZ,4,7),SNQD(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer k
      real a11,a22,d2t,factor,a2200,a2211
	save
c

      if(ijkadv.eq.1) d2t=d2t/2.

       dt_stat=ndt_stat

	factor=dt_stat/d2t
	print *,'(1) dt_stat,d2t,factor=',dt_stat,d2t,factor

       a22=1./float(nx-2)
       a11=dt_stat
       a2200=a22*factor
       a2211=a22*a11
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.003e7

      do k=2,nz-1
       pir=1./PI(K)
       avc=(alv/cp)*pir
       afc=(alf/cp)*pir
       asc=(als/cp)*pir

        acoc(k)=acoc(k)+coc(k,1,4)*a2200
        acoe(k)=acoe(k)+coe(k,1,4)*a2200
        agaut(k)=agaut(k)+(gaut(k,1,4)+tdw(k,1,4))*a2200
        asdep(k)=asdep(k)+sdep(k,1,4)*a2200
        agsub(k)=agsub(k)+gsub(k,1,4)*a2200

        apern(k)=apern(k)+pern(k,1,4)*a2200
        agwet(k)=agwet(k)+gwet(k,1,4)*a2200
        afgr(k)=afgr(k)+gfr(k,1,4)*a2200

	amelt(k)=amelt(k)+(smlt(k,1,4)+gmlt(k,1,4)+ssub(k,1,4))*a2200

        aother(k)=aother(k)+(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     1                      tmlt(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     2                      sacr(k,1,4)-d2rs(k,1,4)+thom(k,1,4))*a2211

        total_mq(k)=total_mq(k)+(coc(k,1,4)-coe(k,1,4)
     1            +gaut(k,1,4)+tdw(k,1,4)+sdep(k,1,4)+gsub(k,1,4)
     2            -pern(k,1,4)
     2            -gwet(k,1,4)-gfr(k,1,4))*a2200       


        total_mt(k)=total_mt(k)+factor*(avc*a22*(coc(k,1,4)-coe(k,1,4))
     1              +asc*a22*(gaut(k,1,4)+tdw(k,1,4)+sdep(k,1,4)
     2                       +gsub(k,1,4)
     2                    -pern(k,1,4)-gwet(k,1,4)-gfr(k,1,4))
     3              -afc*a22*(smlt(k,1,4)+gmlt(k,1,4)+ssub(k,1,4))) 
     4              +afc*a22*(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     5                    tmlt(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     6                    sacr(k,1,4)-d2rs(k,1,4)+thom(k,1,4))*a11

        avelw(k)=avelw(k)+srlw(k,1,4)*a2211
        avesw(k)=avesw(k)+srsw(k,1,4)*a2211
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end do

      return
      end


      subroutine budget_mt_rh (ijkadv,ndt_stat,d2t,dt_stat,
     1                          total_mq,total_mt)
c     modified for TOGA COARE and GATE
      parameter (NX=514,NZ=43)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension total_mq(nz),total_mt(nz)


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

      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNQH(NZ,4,7),
     1 SNQV(NZ,4,7),SNQD(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer k
	real d2t,factor
	save

      if(ijkadv.eq.1)  d2t=d2t/2.

       dt_stat=ndt_stat

        factor=dt_stat/d2t
        print *,'(2) dt_stat,d2t,factor=',dt_stat,d2t,factor

       a22=1./(nx-2)
       a11=dt_stat
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.003e7
       avc=alv/cp
       afc=alf/cp
       asc=als/cp
      do k=2,nz-1

        total_mq(k)=(coc(k,1,4)-coe(k,1,4)
     1            +gaut(k,1,4)+sdep(k,1,4)+gsub(k,1,4)+tdw(k,1,4)
     2            -pern(k,1,4)
     2            -gwet(k,1,4)-gfr(k,1,4))*a22*factor 
c
        total_mt(k)=factor*(avc*a22*(coc(k,1,4)-coe(k,1,4)) 
     1              +asc*a22*(gaut(k,1,4)+tdw(k,1,4)+sdep(k,1,4)
     2                       +gsub(k,1,4)
     2                    -pern(k,1,4)-gwet(k,1,4)-gfr(k,1,4))
     3              -afc*a22*(smlt(k,1,4)+gmlt(k,1,4)+ssub(k,1,4))) 
     4              +afc*a22*(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     5                    tmlt(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     6                    sacr(k,1,4)-d2rs(k,1,4)+thom(k,1,4))*a11


      end do
      return
      end

      subroutine budget_m_lin (ijkadv,ndt_stat,d2t)
      parameter (NX=514,NZ=43)
c     modified for TOGA COARE and GATE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/mbudget/ acoc(nz),acoe(nz),agaut(nz),asdep(nz),agsub(nz),
     1   apern(nz),aqls(nz),altrans(nz),actrans(nz),total_mq(nz),
     2   agwet(nz),afgr(nz),avelw(nz),avesw(nz),amelt(nz),aother(nz),
     3   total_mt(nz)
      common/tqave/ tavet(nz),qavet(nz)

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

      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNQH(NZ,4,7),
     1 SNQV(NZ,4,7),SNQD(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common/b6/ tls(nz),qls(nz),fd(nz),fe(nz),p0(nz),pi(nz),f0(nz),
     1   st(nz),sv(nz),sq(nz),sc(nz),se(nz),sqa(nz),ri(nx),ar(nx),rx(nx)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer k
      real a11,a22,d2t,factor,a2200,a2211
	save
c

      if(ijkadv.eq.1) d2t=d2t/2.

       dt_stat=ndt_stat

	factor=dt_stat/d2t
	print *,'(1) dt_stat,d2t,factor=',dt_stat,d2t,factor

       a22=1./float(nx-2)
       a11=dt_stat
       a2200=a22*factor
       a2211=a22*a11
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.003e7

      do k=2,nz-1
       pir=1./PI(K)
       avc=(alv/cp)*pir
       afc=(alf/cp)*pir
       asc=(als/cp)*pir

        acoc(k)=acoc(k)+coc(k,1,4)*a2200
        acoe(k)=acoe(k)+(coe(k,1,4)+erns(k,1,4))*a2200
        agaut(k)=agaut(k)+d2rs(k,1,4)*a2200
        asdep(k)=asdep(k)+sdep(k,1,4)*a2200
        agsub(k)=agsub(k)+thom(k,1,4)*a2200

        apern(k)=apern(k)+pern(k,1,4)*a2200
        agwet(k)=agwet(k)+gsub(k,1,4)*a2200
        afgr(k)=afgr(k)-ssub(k,1,4)*a2200

	amelt(k)=amelt(k)+(smlt(k,1,4)+gmlt(k,1,4)+tmlt(k,1,4))*a2200

        aother(k)=aother(k)+(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     1                      tdw(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     2                      sacr(k,1,4)+gfr(k,1,4)+wgrs(k,1,4))*a2211

        total_mq(k)=total_mq(k)+(coc(k,1,4)-coe(k,1,4)-erns(k,1,4)
     1            +d2rs(k,1,4)+thom(k,1,4)+sdep(k,1,4)+ssub(k,1,4)
     2            -pern(k,1,4)-gsub(k,1,4))*a2200       


        total_mt(k)=total_mt(k)+factor*(avc*a22*(coc(k,1,4)-coe(k,1,4)
     1                                          -erns(k,1,4))
     1              +asc*a22*(d2rs(k,1,4)+thom(k,1,4)+sdep(k,1,4)
     2                       +ssub(k,1,4)
     2                    -pern(k,1,4)-gsub(k,1,4))
     3              -afc*a22*(smlt(k,1,4)+gmlt(k,1,4)+tmlt(k,1,4))) 
     4              +afc*a22*(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     5                    tdw(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     6                    sacr(k,1,4)+gfr(k,1,4)+wgrs(k,1,4))*a11

        avelw(k)=avelw(k)+srlw(k,1,4)*a2211
        avesw(k)=avesw(k)+srsw(k,1,4)*a2211
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end do

      return
      end


      subroutine budget_mt_lin (ijkadv,ndt_stat,d2t,dt_stat,
     1                          total_mq,total_mt)
c     modified for TOGA COARE and GATE
      parameter (NX=514,NZ=43)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension total_mq(nz),total_mt(nz)


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

      COMMON/BSTS2/ SNTH(NZ,4,7),SNTV(NZ,4,7),SNTD(NZ,4,7),SNQH(NZ,4,7),
     1 SNQV(NZ,4,7),SNQD(NZ,4,7),SNHH(NZ,4,7),SNHV(NZ,4,7),SNHD(NZ,4,7)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer k
	real d2t,factor
	save

      if(ijkadv.eq.1)  d2t=d2t/2.

       dt_stat=ndt_stat

        factor=dt_stat/d2t
        print *,'(2) dt_stat,d2t,factor=',dt_stat,d2t,factor

       a22=1./(nx-2)
       a11=dt_stat
       alv=2.5e10
       alf=3.336e9
       als=2.8336e10
       cp=1.003e7
       avc=alv/cp
       afc=alf/cp
       asc=als/cp
      do k=2,nz-1

        total_mq(k)=(coc(k,1,4)-coe(k,1,4)-erns(k,1,4)
     1            +thom(k,1,4)+sdep(k,1,4)+d2rs(k,1,4)+ssub(k,1,4)
     2            -pern(k,1,4)-gsub(k,1,4))*a22*factor 
c
        total_mt(k)=factor*(avc*a22*(coc(k,1,4)-coe(k,1,4)-erns(k,1,4)) 
     1              +asc*a22*(thom(k,1,4)+d2rs(k,1,4)+sdep(k,1,4)
     2                       +ssub(k,1,4)
     2                    -pern(k,1,4)-gsub(k,1,4))
     3              -afc*a22*(smlt(k,1,4)+gmlt(k,1,4)+tmlt(k,1,4))) 
     4              +afc*a22*(sacw(k,1,4)+tacr(k,1,4)+sfw(k,1,4)+
     5                    tdw(k,1,4)+gacw(k,1,4)+gacr(k,1,4)+
     6                    sacr(k,1,4)+gfr(k,1,4)+wgrs(k,1,4))*a11


      end do
      return
      end
C
CC
C
      SUBROUTINE RINIT_UV
      PARAMETER (NX=514,NZ=43)
      PARAMETER (NZ94=4*9*NZ,NZ24=2*4*NZ)
      real            V1(NX,NZ),V2(NX,NZ),V3(NX,NZ),
     1                V4(NX,NZ),V5(NX,NZ),V6(NX,NZ),
     2                V7(NX,NZ),V8(NX,NZ)
      COMMON/DEBUG_V/ V1,V2,V3,V4,V5,V6,V7,V8

      real            U1(NX,NZ),U2(NX,NZ),U3(NX,NZ),
     1                U4(NX,NZ),U5(NX,NZ),U6(NX,NZ),
     2                U7(NX,NZ),U8(NX,NZ)
      COMMON/DEBUG_U/ U1,U2,U3,U4,U5,U6,U7,U8

      real             SDIFF_VX(NZ94)
      COMMON/DEBUG_SV/ SDIFF_VX

      real             SDIFF_UX(NZ94)
      COMMON/DEBUG_SU/ SDIFF_UX

      real              SSDIFF_VX(NZ94)
      COMMON/DEBUG_SVS/ SSDIFF_VX

      real              SSDIFF_UX(NZ94)
      COMMON/DEBUG_SUS/ SSDIFF_UX

      real        UB_MEAN(NZ24)
      COMMON/SUV/ UB_MEAN

      real U9(NX,NZ),V9(NX,NZ)
      COMMON/DEBUG_UUVV/ U9,V9

      real        SUB_MEAN(NZ24)
      COMMON/SSUV/ SUB_MEAN

      DO K=1,NZ
         DO I=1,NX
            V1(I,K)=0.
            V2(I,K)=0.
            V3(I,K)=0.
            V4(I,K)=0.
            V5(I,K)=0.
            V6(I,K)=0.
            V7(I,K)=0.
            V8(I,K)=0.
C
            U1(I,K)=0.
            U2(I,K)=0.
            U3(I,K)=0.
            U4(I,K)=0.
            U5(I,K)=0.
            U6(I,K)=0.
            U7(I,K)=0.
            U8(I,K)=0.

            U9(I,K)=0.
            V9(I,K)=0.
         ENDDO
      ENDDO
      DO K=1,NZ94
         SDIFF_VX(K)=0.
         SDIFF_UX(K)=0.
         SSDIFF_VX(K)=0.
         SSDIFF_UX(K)=0.
      ENDDO
      DO K=1,NZ24
         UB_MEAN(K)=0.
         SUB_MEAN(K)=0.
      ENDDO
C
      RETURN
      END
c
      SUBROUTINE RINIT_UV2
      PARAMETER (NX=514,NZ=43)
      PARAMETER (NZ94=4*9*NZ,NZ24=2*4*NZ)
      real            V1(NX,NZ),V2(NX,NZ),V3(NX,NZ),
     1                V4(NX,NZ),V5(NX,NZ),V6(NX,NZ),
     2                V7(NX,NZ),V8(NX,NZ)
      COMMON/DEBUG_V/ V1,V2,V3,V4,V5,V6,V7,V8

      real            U1(NX,NZ),U2(NX,NZ),U3(NX,NZ),
     1                U4(NX,NZ),U5(NX,NZ),U6(NX,NZ),
     2                U7(NX,NZ),U8(NX,NZ)
      COMMON/DEBUG_U/ U1,U2,U3,U4,U5,U6,U7,U8

      real U9(NX,NZ),V9(NX,NZ)
      COMMON/DEBUG_UUVV/ U9,V9

      DO K=1,NZ
         DO I=1,NX
            V1(I,K)=0.
            V2(I,K)=0.
            V3(I,K)=0.
            V4(I,K)=0.
            V5(I,K)=0.
            V6(I,K)=0.
            V7(I,K)=0.
            V8(I,K)=0.
C
            U1(I,K)=0.
            U2(I,K)=0.
            U3(I,K)=0.
            U4(I,K)=0.
            U5(I,K)=0.
            U6(I,K)=0.
            U7(I,K)=0.
            U8(I,K)=0.

            U9(I,K)=0.
            V9(I,K)=0.
         ENDDO
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE WTAP_UV
      PARAMETER (NX=514,NZ=43)
      real            DIFF_VX(NX,NZ),DIFF_VZ(NX,NZ),GRID_VX(NX,NZ),
     1                GRID_VZ(NX,NZ),DIFF_NV(NX,NZ),V_LARGE(NX,NZ),
     2                PRE_V(NX,NZ),DT_VWIND(NX,NZ)
      COMMON/DEBUG_V/ DIFF_VX,DIFF_VZ,GRID_VX,
     1                GRID_VZ,DIFF_NV,V_LARGE,
     2                PRE_V,DT_VWIND

      real            DIFF_UX(NX,NZ),DIFF_UZ(NX,NZ),GRID_UX(NX,NZ),
     1                GRID_UZ(NX,NZ),DIFF_NU(NX,NZ),U_LARGE(NX,NZ),
     2                PRE_U(NX,NZ),DT_UWIND(NX,NZ)
      COMMON/DEBUG_U/ DIFF_UX,DIFF_UZ,GRID_UX,
     1                GRID_UZ,DIFF_NU,U_LARGE,
     2                PRE_U,DT_UWIND

C
      REAL BNXNZ(NX,NZ)
C
       CALL WTAP2 ( DIFF_VX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_VZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_VX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_VZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_NV,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( V_LARGE,BNXNZ,NX,NZ,4)
       CALL WTAP2 (   PRE_V,BNXNZ,NX,NZ,4)
       CALL WTAP2 (DT_VWIND,BNXNZ,NX,NZ,4)
CCCCCCC
       CALL WTAP2 ( DIFF_UX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_UZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_UX,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( GRID_UZ,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( DIFF_NU,BNXNZ,NX,NZ,4)
       CALL WTAP2 ( U_LARGE,BNXNZ,NX,NZ,4)
       CALL WTAP2 (   PRE_U,BNXNZ,NX,NZ,4)
       CALL WTAP2 (DT_UWIND,BNXNZ,NX,NZ,4)
      RETURN
      END
C
      SUBROUTINE WTAP_UV_S
      PARAMETER (NZ=43)
      real             SDIFF_VX(NZ,4),SDIFF_VZ(NZ,4),SGRID_VX(NZ,4),
     1                 SGRID_VZ(NZ,4),SDIFF_NV(NZ,4),SV_LARGE(NZ,4),
     2                 SPRE_V(NZ,4),SDT_V(NZ,4),TN_VWIND(NZ,4)
      COMMON/DEBUG_SV/ SDIFF_VX,SDIFF_VZ,SGRID_VX,
     1                 SGRID_VZ,SDIFF_NV,SV_LARGE,
     2                 SPRE_V,SDT_V,TN_VWIND

      real             SDIFF_UX(NZ,4),SDIFF_UZ(NZ,4),SGRID_UX(NZ,4),
     1                 SGRID_UZ(NZ,4),SDIFF_NU(NZ,4),SU_LARGE(NZ,4),
     2                 SPRE_U(NZ,4),SDT_U(NZ,4),TN_UWIND(NZ,4)
      COMMON/DEBUG_SU/ SDIFF_UX,SDIFF_UZ,SGRID_UX,
     1                 SGRID_UZ,SDIFF_NU,SU_LARGE,
     2                 SPRE_U,SDT_U,TN_UWIND

      real        UB_MEAN(NZ,4),VB_MEAN(NZ,4)
      COMMON/SUV/ UB_MEAN,VB_MEAN

      real             SSDIFF_VX(NZ,4),SSDIFF_VZ(NZ,4),SSGRID_VX(NZ,4),
     1                 SSGRID_VZ(NZ,4),SSDIFF_NV(NZ,4),SSV_LARGE(NZ,4),
     2                 SSPRE_V(NZ,4),SSDT_V(NZ,4),STN_VWIND(NZ,4)
      COMMON/DEBUG_SVS/ SSDIFF_VX,SSDIFF_VZ,SSGRID_VX,
     1                  SSGRID_VZ,SSDIFF_NV,SSV_LARGE,
     2                  SSPRE_V,SSDT_V,STN_VWIND

      real             SSDIFF_UX(NZ,4),SSDIFF_UZ(NZ,4),SSGRID_UX(NZ,4),
     1                 SSGRID_UZ(NZ,4),SSDIFF_NU(NZ,4),SSU_LARGE(NZ,4),
     2                 SSPRE_U(NZ,4),SSDT_U(NZ,4),STN_UWIND(NZ,4)
      COMMON/DEBUG_SUS/ SSDIFF_UX,SSDIFF_UZ,SSGRID_UX,
     1                  SSGRID_UZ,SSDIFF_NU,SSU_LARGE,
     2                  SSPRE_U,SSDT_U,STN_UWIND

      real        SUB_MEAN(NZ,4),SVB_MEAN(NZ,4)
      COMMON/SSUV/ SUB_MEAN,SVB_MEAN

C
      REAL BNZ4(NZ,4)
C
      CALL WTAP2 (SDIFF_VX,BNZ4,NZ,4,4)
      CALL WTAP2 (SDIFF_VZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SGRID_VX,BNZ4,NZ,4,4)
      CALL WTAP2 (SGRID_VZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SDIFF_NV,BNZ4,NZ,4,4)
      CALL WTAP2 (SV_LARGE,BNZ4,NZ,4,4)
      CALL WTAP2 (  SPRE_V,BNZ4,NZ,4,4)
      CALL WTAP2 (   SDT_V,BNZ4,NZ,4,4)
      CALL WTAP2 (TN_VWIND,BNZ4,NZ,4,4)
CCCCCC
      CALL WTAP2 (SDIFF_UX,BNZ4,NZ,4,4)
      CALL WTAP2 (SDIFF_UZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SGRID_UX,BNZ4,NZ,4,4)
      CALL WTAP2 (SGRID_UZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SDIFF_NU,BNZ4,NZ,4,4)
      CALL WTAP2 (SU_LARGE,BNZ4,NZ,4,4)
      CALL WTAP2 (  SPRE_U,BNZ4,NZ,4,4)
      CALL WTAP2 (   SDT_U,BNZ4,NZ,4,4)
      CALL WTAP2 (TN_UWIND,BNZ4,NZ,4,4)
cccccc
      CALL WTAP2 ( UB_MEAN,BNZ4,NZ,4,4)
      CALL WTAP2 ( VB_MEAN,BNZ4,NZ,4,4)
C
CC
C
      CALL WTAP2 (SSDIFF_VX,BNZ4,NZ,4,4)
      CALL WTAP2 (SSDIFF_VZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SSGRID_VX,BNZ4,NZ,4,4)
      CALL WTAP2 (SSGRID_VZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SSDIFF_NV,BNZ4,NZ,4,4)
      CALL WTAP2 (SSV_LARGE,BNZ4,NZ,4,4)
      CALL WTAP2 (  SSPRE_V,BNZ4,NZ,4,4)
      CALL WTAP2 (   SSDT_V,BNZ4,NZ,4,4)
      CALL WTAP2 (STN_VWIND,BNZ4,NZ,4,4)
CCCCCC
      CALL WTAP2 (SSDIFF_UX,BNZ4,NZ,4,4)
      CALL WTAP2 (SSDIFF_UZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SSGRID_UX,BNZ4,NZ,4,4)
      CALL WTAP2 (SSGRID_UZ,BNZ4,NZ,4,4)
      CALL WTAP2 (SSDIFF_NU,BNZ4,NZ,4,4)
      CALL WTAP2 (SSU_LARGE,BNZ4,NZ,4,4)
      CALL WTAP2 (  SSPRE_U,BNZ4,NZ,4,4)
      CALL WTAP2 (   SSDT_U,BNZ4,NZ,4,4)
      CALL WTAP2 (STN_UWIND,BNZ4,NZ,4,4)
cccccc
      CALL WTAP2 ( SUB_MEAN,BNZ4,NZ,4,4)
      CALL WTAP2 ( SVB_MEAN,BNZ4,NZ,4,4)
cccccc
C
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wstrad2 (itape)  ! newly added! 8/21/01 shie
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   8/21/01 by shie, write more radiation infromation to restart file

      parameter (NX=514,lay=88)

      common/srflx/ sfir(lay,4),sfsw(lay,4),shir(lay,4),shsw(lay,4),
     $        salpha(4),si(4)

      common/surface/ sun_4(nx,4)

      save
      if(itape.gt.6) go to 100
       write(itape) sfir
       write(itape) sfsw
       write(itape) shir
       write(itape) shsw
       write(itape) salpha
       write(itape) si
       write(itape) sun_4
      return
  100  read(itape) sfir
       read(itape) sfsw
       read(itape) shir
       read(itape) shsw
       read(itape) salpha
       read(itape) si
       read(itape) sun_4
C
      RETURN
      END

