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
