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
