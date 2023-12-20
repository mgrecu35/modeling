
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
