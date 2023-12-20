

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
