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
