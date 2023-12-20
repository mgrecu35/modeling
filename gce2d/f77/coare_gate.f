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
      real rlh(nx,nz)
      common /lhblock/ rlh
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
      real dtint
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
      !print*, dz0
      !stop
      OPEN( 4,FILE='data.3min0',FORM='UNFORMATTED')
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
      dtint=300.
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
      print*, dz
      print*, 'dz'
      stop
      nran=0
      rm=rm0/150.
      do 107 i=2,iles
 107    rx(i)= (ran3(ised1)-0.5)*rm
      nran=nran+1
c     ***********************
      if (irs.eq.0 .and. n.eq.1) go to 6666
      if (irs.eq.0 .and. n.eq.0) go to 5555
 
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
      call initlh()
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



C     TAO 12-10-97
        if(icps_bm .eq. 1) then
        endif
C     TAO 12-10-97


cc    *********************************
 1111 if(isec.lt.iwtp) go to 2222



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

c     ************************
 9999 continue
      if(isec/int(dtint)*int(dtint).eq.isec) then
         write(*,*) 'LHaverage'
        call average(dtint)
      endif
cwrite(*,*) 'ISEC=',isec
      if(isec/180*180.eq.isec) then	
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
       call wtap2 (rlh,bnxnz,nx,nz,4)
       call wtap1 (rho,bnz,nz,4)
       call wtap1 (ta,bnz,nz,4)
       call wtap1 (qa,bnz,nz,4)
       call wtap1 (pi,bnz,nz,4)
       call wtap1 (p0,bnz,nz,4)
      endif 
      if(isec/int(dtint)*int(dtint).eq.isec) then
          call initlh()
      endif

 9990 if(isec.lt.iend) go to 3333
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

      if(icps_bm .eq. 1) then

      endif
     
      if(istaton.eq.1)then

      endif

C     Write to restart file

        rdt=dt

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
